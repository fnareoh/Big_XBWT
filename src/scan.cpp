#include "parameters.hpp"
#include <algorithm>
#include <assert.h>
#include <ctime>
#include <errno.h>
#include <fstream>
#include <iomanip>
#include <map>
#include <random>
#include <sstream>
#include <stdexcept>
#include <stdint.h>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>
extern "C" {
#include "../external/utils.h"
}
#ifdef BAM_READER
// To read the bam format
#include "api/BamReader.h"
#include "api/BamWriter.h"
using namespace BamTools;
#endif

using namespace std;
// using namespace __gnu_cxx;

// =============== algorithm limits ===================
// maximum number of distinct words
#define MAX_DISTINCT_WORDS (INT32_MAX - 1)
typedef uint32_t word_int_t;
// maximum number of occurrences of a single word
#define MAX_WORD_OCC (UINT32_MAX)
typedef uint32_t occ_int_t;

// values of the wordFreq map: word, its number of occurrences, and its rank
struct word_stats {
  string str;
  occ_int_t occ;
  word_int_t rank = 0;
};

// -----------------------------------------------------------------
// class to maintain a window in a string and its KR fingerprint
struct KR_window {
  int wsize;
  int *window;
  int asize;
  const uint64_t prime = 1999999973;
  uint64_t hash;
  uint64_t tot_char;
  uint64_t asize_pot; // asize^(wsize-1) mod prime

  KR_window(int w) : wsize(w) {
    asize = 256;
    asize_pot = 1;
    for (int i = 1; i < wsize; i++)
      asize_pot =
          (asize_pot * asize) % prime; // ugly linear-time power algorithm
    // alloc and clear window
    window = new int[wsize];
    reset();
  }

  // init window, hash, and tot_char
  void reset() {
    for (int i = 0; i < wsize; i++)
      window[i] = 0;
    // init hash value and related values
    hash = tot_char = 0;
  }

  uint64_t addchar(int c) {
    int k = tot_char++ % wsize;
    // complex expression to avoid negative numbers
    hash += (prime -
             (window[k] * asize_pot) % prime); // remove window[k] contribution
    hash = (asize * hash + c) % prime;         //  add char i
    window[k] = c;
    // cerr << get_window() << " ~~ " << window << " --> " << hash << endl;
    return hash;
  }
  // debug only
  string get_window() {
    string w = "";
    int k = (tot_char - 1) % wsize;
    for (int i = k + 1; i < k + 1 + wsize; i++)
      w.append(1, window[i % wsize]);
    return w;
  }

  ~KR_window() { delete[] window; }
};
// -----------------------------------------------------------

template <typename T = uint32_t> void write_binary(T item, ofstream &stream) {
  stream.write((char *)&item, sizeof(item));
}

static void save_update_word(Args &arg, string &w,
                             map<uint64_t, word_stats> &freq,
                             FILE *tmp_parse_file, bool is_ref,
                             vector<pair<uint64_t, uint64_t>> &start_phrase,
                             FILE *sa, uint64_t &pos);

// compute 64-bit KR hash of a string
// to avoid overflows in 64 bit aritmethic the prime is taken < 2**55
uint64_t kr_hash(string s) {
  uint64_t hash = 0;
  // const uint64_t prime = 3355443229;     // next prime(2**31+2**30+2**27)
  // const uint64_t prime = 27162335252586509; // next prime (2**54 + 2**53 +
  // 2**47 + 2**13)
  for (size_t k = 0; k < s.size(); k++) {
    int c = (unsigned char)s[k];
    assert(c >= 0 && c < 256);
    hash = (256 * hash + c) % PRIME; //  add char k
  }
  return hash;
}

// save current word in the freq map and update it leaving only the
// last minsize chars which is the overlap with next word
static void save_update_word(Args &arg, string &w,
                             map<uint64_t, word_stats> &freq,
                             FILE *tmp_parse_file, bool is_ref,
                             vector<pair<uint64_t, uint64_t>> &start_phrase,
                             FILE *sa, uint64_t &pos) {
  size_t minsize = arg.w;
  if (w.size() <= minsize)
    return;
  cout << "word: " << w << endl;
  // save overlap
  string overlap(w.substr(w.size() - minsize)); // keep last minsize chars

  // get the hash value and write it to the temporary parse file
  uint64_t hash = kr_hash(w);
  if (fwrite(&hash, sizeof(hash), 1, tmp_parse_file) != 1)
    die("parse write error");
  // recall the start of phrases in the reference so we can use them to extend
  // the reads.
  if (is_ref) {
    if (pos == 0) // we don't extend with the dollars
      start_phrase.push_back(make_pair(pos, hash));
    else
      start_phrase.push_back(make_pair(pos, hash));
  }

  // update frequency table for current hash
  if (freq.find(hash) == freq.end()) {
    freq[hash].occ = 1; // new hash
    freq[hash].str = w;
  } else {
    freq[hash].occ += 1; // known hash
    if (freq[hash].occ <= 0) {
      cerr << "Emergency exit! Maximum # of occurence of dictionary word (";
      cerr << MAX_WORD_OCC << ") exceeded\n";
      exit(1);
    }
    if (freq[hash].str != w) {
      cerr << "Emergency exit! Hash collision for strings:\n";
      cerr << freq[hash].str << "\n  vs\n" << w << endl;
      exit(1);
    }
  }
  cout << "nb occ word: " << freq[hash].occ << endl;

  // update sa files
  // compute ending position +1 of current word and write it to sa file
  // pos is the ending position+1 of the previous word and is updated here
  if (pos == 0)
    // minus minsize because of the w initials $ of the first word
    pos = w.size() - minsize;
  else
    pos += w.size() - minsize;
  if (sa)
    if (fwrite(&pos, IBYTES, 1, sa) != 1)
      die("Error writing to sa info file");

  // keep only the overlapping part of the window
  w.assign(overlap);
}

// prefix free parse of the reference and the reads. w is the window size,
// p is the modulus use a KR-hash as the word ID that is immediately written
// to the temporary parse file : input.parse_old
uint64_t process_file(Args &arg, map<uint64_t, word_stats> &wordFreq) {

  // ----------------------- Parsing the reference  -------------------------o--

  // open a input file for the reference
  string fnam = arg.inputFileName;
  bool is_fasta = false;
  if (fnam.substr(fnam.find_last_of(".") + 1) == "fasta")
    is_fasta = true;
  ifstream f(fnam);

  if (!f.rdbuf()->is_open()) { // is_open does not work on igzstreams
    perror(__func__);
    throw new std::runtime_error("Cannot open input file " + fnam);
  }

  // open the 1st pass parsing file
  FILE *tmp_parse_file = open_aux_file(arg.inputFileName.c_str(), EXTPARS0, "wb");
  // (Starting position of a phrase, hash of this phrase)
  vector<pair<uint64_t, uint64_t>> start_phrase;
  FILE *sa_file = NULL;
  auto limit_file = ofstream(arg.inputFileName + "." + EXTLIM);
  // if requested open file containing the ending position+1 of each word
  if (arg.SAinfo)
    sa_file = open_aux_file(arg.inputFileName.c_str(), EXTSAI, "wb");

  // main loop on the chars of the input file
  cout << "Parsing the genome" << endl;
  char c;
  string line = "";   // storing of a string used for the fasta format
  int index_line = 0; // storing the position in the line
  uint64_t pos = 0; // ending position +1 of previous word in the original text,
                    // used for computing sa_info
  assert(IBYTES <=
         sizeof(pos)); // IBYTES bytes of pos are written to the sa info file
  // init first word in the parsing with w Dollar char (s we don't forget any
  // valuable char at the last step
  string word("");
  word.append(arg.w, Dollar);
  // init empty KR window: constructor only needs window size
  KR_window krw(arg.w);
#ifdef OUTPUT_EXTENDED_READ
  ofstream extended_file;
  extended_file.open(arg.inputFileName + ".extended_input");
  extended_file << arg.w << endl;
#endif
  while (true) {
    if (is_fasta) {
      if ((long)index_line < (long)line.size()) {
        c = line[index_line];
        index_line++;
        if (c == '\n')
          continue;
      } else {
        if (f.eof())
          break;
        getline(f, line);
        index_line = 0;
        if (line[0] == '>') {
          index_line = line.size();
        }
        continue;
      }
    }
    if (!is_fasta && ((c = f.get()) == EOF))
      break;
    if (c == '\n')
      continue;
    word.append(1, c);
#ifdef OUTPUT_EXTENDED_READ
    extended_file << (char)c;
#endif
    uint64_t hash = krw.addchar(c);
    //cout << c << " " << hash << endl;
    if (hash % arg.p == 0) {
      // end of word, save it and write its full hash to the output file
      // cerr << "~"<< c << "~ " << hash << " ~~ <" << word << "> ~~ <" <<
      // krw.get_window() << ">" <<  endl;
      uint32_t l_p_start = 0;
      uint32_t l_p_end = word.size();
      cout << "Ref l_start l_end: " << l_p_start << " " << l_p_end << endl;
      if (word.size() >  (unsigned) arg.w){
        write_binary(l_p_start, limit_file);
        write_binary(l_p_end, limit_file);
      }
      save_update_word(arg, word, wordFreq, tmp_parse_file, true, start_phrase, sa_file,
                       pos);
    }
  }
  // add the last word in the dict
  uint32_t l_p_start = 0;
  uint32_t l_p_end = word.size();
  cout << "Last ref l_start l_end: " << l_p_start << " " << l_p_end << endl;
  if (word.size() > (unsigned) arg.w){
    write_binary(l_p_start, limit_file);
    write_binary(l_p_end, limit_file);
  }
  save_update_word(arg, word, wordFreq, tmp_parse_file, true, start_phrase, sa_file, pos);

  assert(pos == krw.tot_char);
  cout << "Length of the parsed reference: " << start_phrase.size() << endl;
  uint64_t totChar = krw.tot_char;

#ifdef OUTPUT_EXTENDED_READ
  extended_file << endl;
#endif
  // ----------------------- Parsing the reads  -------------------------------
  cout << "Parsing the reads" << endl;

  // open read file
  fnam = arg.ReadsFileName;

  // Boolean to know wheter we are using the bam format
  bool is_bam = false;
  if (fnam.substr(fnam.find_last_of(".") + 1) == "bam")
    is_bam = true;

  // file stream if we are not using the bam format
  ifstream f_r;
  if (!is_bam) {
    f_r.open(fnam);
    if (!f_r.rdbuf()->is_open()) {
      perror(__func__);
      throw new std::runtime_error("Cannot open input file " + fnam);
    }
  }

#ifdef BAM_READER
  // if we are using the bam format get the BamReader
  BamReader reader;
  if (is_bam) {
    if (!reader.Open(fnam)) {
      cerr << "Could not open input BAM file." << endl;
      return 1;
    }
  }
  BamAlignment al;
#endif

  uint64_t pos_read;
  string read;
  // Process each read
  while (true) {
    // ----------------------- Get the input ----------------------------------

    // get input if we are not using the bam format
    if (!is_bam) {
      if (f_r.eof())
        break;
      f_r >> pos_read >> read;
      // if (f_r.eof())
      // break;

    }
#ifdef BAM_READER
    // get input if we are using the bam format
    else {
      if (!reader.GetNextAlignment(al))
        break;
      pos_read = al.Position;
      read = al.QueryBases;
      if (read == "" || read == "*")
        continue;
      if (pos_read == -1)
        continue;
    }
#endif

    // ------------------ Retrieve the extended read --------------------------

    // starting position of the extended read in the parse
    uint64_t r_s_p =
        upper_bound(start_phrase.begin(), start_phrase.end(),
                    make_pair(pos_read, numeric_limits<uint64_t>::max())) -
        start_phrase.begin() - 1;
    // ending position of the extended read in the parse
    uint64_t r_e_p = upper_bound(start_phrase.begin(), start_phrase.end(),
                                 make_pair(pos_read + read.size()-1,
                                           numeric_limits<uint64_t>::max())) -
                     start_phrase.begin() - 1;
    /*cout << "pos_read " << pos_read << endl;
    cout << "pos_read_end " << pos_read + read.size() << endl;
    cout << "r_s_p " << r_s_p << endl;
    cout << "r_e_p " << r_e_p << endl;*/

    assert(r_s_p < start_phrase.size());
    assert(r_e_p < start_phrase.size());
    // phrase we are going to extend the read with, at the front and at the end
    string front_phrase = wordFreq[start_phrase[r_s_p].second].str;
    string back_phrase = wordFreq[start_phrase[r_e_p].second].str;
    /*cout << "start_phrases: " << start_phrase[r_s_p].first << " end_phrase "
         << start_phrase[r_e_p].first + size(back_phrase) << endl;
    cout << "start_phrases r_e_p+1: " << start_phrase[r_e_p + 1].first << endl;
    cout << "len_phrases: " << size(front_phrase) << " " << size(back_phrase)
         << endl;*/

    // The extended read that will be parsed in to phrases
    string read_extanded;
    uint32_t l_start;
    uint32_t l_end;
    if (pos_read + read.size() >
        start_phrase[r_e_p].first + size(back_phrase)) {
      l_start = pos_read + arg.w - start_phrase[r_s_p].first;
      read_extanded = front_phrase.substr(0, l_start) + read;
      l_end = read_extanded.size();
    } else {
      l_start = pos_read + arg.w - start_phrase[r_s_p].first;
      l_end = pos_read + read.size() - start_phrase[r_e_p].first + arg.w;
      // cout << "l_start: " << l_start << endl;
      // cout << "l_end: " << l_end << endl;
      read_extanded =
          front_phrase.substr(0, l_start) + read + back_phrase.substr(l_end);
      l_end = read_extanded.size() - back_phrase.size() + l_end;
      // cout << "final l_end: " << l_end << endl;
    }
    // l_start = l_start-arg.w;
    // l_end = l_end-arg.w;
    cout << "read: " << read << endl;
    cout << "read_extanded: " << read_extanded << endl;
    cout << "l_start: " << l_start << endl;
    cout << "l_end: " << l_end << endl;
#ifdef OUTPUT_EXTENDED_READ
    extended_file << start_phrase[r_s_p].first << " "
                  << read_extanded.substr(10) << endl;
#endif

    // ------------------- Parse the extended read ----------------------------

    // Put a sepatator in the parsing
    uint64_t separator = PRIME + 1;
    if (fwrite(&separator, sizeof(separator), 1, tmp_parse_file) != 1)
      die("parse write error");
    // Write the start of the extended read in the parse
    if (fwrite(&r_s_p, sizeof(r_s_p), 1, tmp_parse_file) != 1)
      die("parse write error");

    // Write the parsing of the extended read
    // init empty KR window: constructor only needs window size
    krw.reset();
    uint64_t i = 0;
    bool after_start = false;
    pos = start_phrase[r_s_p].first;

    word = read_extanded.substr(0, arg.w);
    while ((long)i < (long)arg.w) {
      if (pos != 0) krw.addchar(read_extanded[i]);
      //cout << read_extanded[i] << endl;
      i++;
    }

    while (i < read_extanded.size()) {
      // cout << "i: " << i << " word.size(): " << word.size() << endl;
      word.append(1, read_extanded[i]);
      uint64_t hash = krw.addchar(read_extanded[i]);
      //cout << read_extanded[i] << " " << hash << endl;
      if (hash % arg.p == 0) {
        // end of word, save it and write its full hash to the output file
        uint32_t l_p_start;
        uint32_t l_p_end;
        if (i < l_start)
          l_p_start = word.size();
        else {
          if (!after_start) {
            after_start = true;
            cout << "i: " << i << endl;
            cout << "l_start: " << l_start << endl;
            cout << "word.size(): " << word.size() << endl;
            cout << "pos_read: " << pos_read << endl;
            l_p_start = l_start - i + word.size()-1;
          } else
            l_p_start = 0;
        }
        if (i < l_end)
          l_p_end = word.size();
        else
          l_p_end = word.size() - i + l_end -1;
        /*cout << "i: " << i << endl;
        cout << "word.size(): " << word.size() << endl;*/
        cout << "WRITE l_p_start: " << l_p_start << endl;
        cout << "WRITE l_p_end: " << l_p_end << endl;
        write_binary(l_p_start, limit_file);
        write_binary(l_p_end, limit_file);
        save_update_word(arg, word, wordFreq, tmp_parse_file, false, start_phrase, sa_file,
                         pos);
      }
      i++;
    }
    if (!((start_phrase[r_s_p].first == 0 && pos + 1 == read_extanded.size()) ||
          (pos - start_phrase[r_s_p].first + arg.w == read_extanded.size()))) {
      // If we did not finished on a triggering substring, we add the last
      // phrase.
      uint32_t l_p_start;
      uint32_t l_p_end;
      if (i < l_start)
        l_p_start = word.size();
      else {
        if (!after_start) {
          after_start = true;
          l_p_start = l_start - i + word.size();
        } else
          l_p_start = 0;
      }
      if (i < l_end)
        l_p_end = word.size();
      else
        l_p_end = word.size() - i + l_end;
      /*cout << "i: " << i << endl;
      cout << "word.size(): " << word.size() << endl;*/
      cout << "LAST-WRITE end phrase l_p_start: " << l_p_start << endl;
      cout << "LAST-WRITE end phrase l_p_end: " << l_p_end << endl;
      write_binary(l_p_start, limit_file);
      write_binary(l_p_end, limit_file);
      save_update_word(arg, word, wordFreq, tmp_parse_file, false, start_phrase, sa_file,
                       pos);
    }
    totChar += read_extanded.size() - arg.w;
  }
#ifdef OUTPUT_EXTENDED_READ
  extended_file.close();
#endif
  // close input and output files
  if (sa_file)
    if (fclose(sa_file) != 0)
      die("Error closing SA file");
  if (fclose(tmp_parse_file) != 0)
    die("Error closing parse file");
  f.close();
  return totChar;
}

// function used to compare two string pointers
bool pstringCompare(const string *a, const string *b) { return *a <= *b; }

// given the sorted dictionary and the frequency map write the dictionary and
// occ files also compute the 1-based rank for each hash
void writeDictOcc(Args &arg, map<uint64_t, word_stats> &wfreq,
                  vector<string *> &sortedDict) {
  assert(sortedDict.size() == wfreq.size());
  FILE *fdict, *focc = NULL;
  // open dictionary and occ files
  fdict = open_aux_file(arg.inputFileName.c_str(), EXTDICT, "wb");
  focc = open_aux_file(arg.inputFileName.c_str(), EXTOCC, "wb");

  word_int_t wrank = 1; // current word rank (1 based)
  for (auto x :
       sortedDict) { // *x is the string representing the dictionary word
    char *word = (*x).data();  // current dictionary word
    size_t len = (*x).size();  // offset and length of word
    reverse(word, word + len); // Reverse back to get the hash
    cout << "phrase:" << word << endl;
    assert(len > (size_t)arg.w);
    uint64_t hash = kr_hash(*x);
    auto &wf = wfreq.at(hash);
    assert(wf.occ > 0);
    // Reverse back to write the words in reverse in the dict for the last step
    reverse(word, word + len);
    size_t s = fwrite(word, 1, len, fdict);
    if (s != len)
      die("Error writing to DICT file");
    if (fputc(EndOfWord, fdict) == EOF)
      die("Error writing EndOfWord to DICT file");
    s = fwrite(&wf.occ, sizeof(wf.occ), 1, focc);
    if (s != 1)
      die("Error writing to OCC file");
    assert(wf.rank == 0);
    wf.rank = wrank++;
  }
  if (fputc(EndOfDict, fdict) == EOF)
    die("Error writing EndOfDict to DICT file");
  if (fclose(focc) != 0)
    die("Error closing OCC file");
  if (fclose(fdict) != 0)
    die("Error closing DICT file");
}

void remapParse(Args &arg, map<uint64_t, word_stats> &wfreq) {
  // open parse files. the old parse can be stored in a single file or in
  // multiple files
  mFile *moldp = mopen_aux_file(arg.inputFileName.c_str(), EXTPARS0, arg.th);
  FILE *newp = open_aux_file(arg.inputFileName.c_str(), EXTPARSE, "wb");

  // recompute occ as an extra check
  vector<occ_int_t> occ(wfreq.size() + 1, 0); // ranks are zero based
  uint64_t hash;
  uint64_t separator = PRIME + 1; // separator in parse_old
  uint32_t alphabet_parse = wfreq.size();
  size_t s = fwrite(&alphabet_parse, sizeof(alphabet_parse), 1, newp);
  if (s != 1)
    die("Error writing to new parse file");
  while (true) {
    s = mfread(&hash, sizeof(hash), 1, moldp);
    if (s == 0)
      break;
    if (s != 1)
      die("Unexpected parse EOF");
    if (hash == separator) {
      uint32_t sep_newp = wfreq.size() + 1; // separator in parse
      s = fwrite(&sep_newp, sizeof(sep_newp), 1, newp);
      if (s != 1)
        die("Error writing to new parse file");
      // Read the alignment and rewrite it the same way
      s = mfread(&hash, sizeof(hash), 1, moldp);
      if (fwrite(&hash, sizeof(hash), 1, newp) != 1)
        die("Error writing to new parse file");
    } else {
      // output the remaped character of the parse
      uint32_t rank = wfreq.at(hash).rank;
      occ[rank]++;
      s = fwrite(&rank, sizeof(rank), 1, newp);
      if (s != 1)
        die("Error writing to new parse file");
    }
  }
  if (fclose(newp) != 0)
    die("Error closing new parse file");
  if (mfclose(moldp) != 0)
    die("Error closing old parse segment");
  // check old and recomputed occ coincide
  for (auto &x : wfreq)
    assert(x.second.occ == occ[x.second.rank]);
}

int main(int argc, char **argv) {
  // translate command line parameters
  Args arg;
  scan_parseargs(argc, argv, arg);
  cout << "Windows size: " << arg.w << endl;
  cout << "Stop word modulus: " << arg.p << endl;

  // measure elapsed wall clock time
  time_t start_main = time(NULL);
  time_t start_wc = start_main;
  // init sorted map counting the number of occurrences of each word
  map<uint64_t, word_stats> wordFreq;
  uint64_t totChar;

  // ---------------------- proscessing input files ---------------------------
  try {
    if (arg.th == 0)
      totChar = process_file(arg, wordFreq);
    else {
      cerr << "Sorry, this is the no-threads executable and you requested "
           << arg.th << " threads\n";
      exit(EXIT_FAILURE);
    }
  } catch (const std::bad_alloc &) {
    cout << "Out of memory (parsing phase)... emergency exit\n";
    die("bad alloc exception");
  }
  // first report
  uint64_t totDWord = wordFreq.size();
  cout << "Total input symbols: " << totChar << endl;
  cout << "Found " << totDWord << " distinct words" << endl;
  cout << "Parsing took: " << difftime(time(NULL), start_wc)
       << " wall clock seconds\n";
  // check # distinct words
  if (totDWord > MAX_DISTINCT_WORDS) {
    cerr << "Emergency exit! The number of distinc words (" << totDWord
         << ")\n";
    cerr << "is larger than the current limit (" << MAX_DISTINCT_WORDS << ")\n";
    exit(1);
  }

  // -------------- second pass
  start_wc = time(NULL);
  // create array of dictionary words
  vector<string *> dictArray;
  dictArray.reserve(totDWord);
  // fill array
  uint64_t sumLen = 0;
  uint64_t totWord = 0;
  for (auto &x : wordFreq) {
    sumLen += x.second.str.size();
    totWord += x.second.occ;
    // reverse so that we sort in colexicographic order
    reverse(x.second.str.begin(), x.second.str.end());
    dictArray.push_back(&x.second.str);
  }
  assert(dictArray.size() == totDWord);
  cout << "Sum of lenghts of dictionary words: " << sumLen << endl;
  cout << "Total number of words: " << totWord << endl;
  // sort dictionary
  sort(dictArray.begin(), dictArray.end(), pstringCompare);
  // write plain dictionary and occ file, also compute rank for each hash
  cout << "Writing plain dictionary and occ file\n";
  writeDictOcc(arg, wordFreq, dictArray);
  dictArray.clear(); // reclaim memory
  cout << "Dictionary construction took: " << difftime(time(NULL), start_wc)
       << " wall clock seconds\n";

  // remap parse file
  start_wc = time(NULL);
  cout << "Generating remapped parse file\n";
  remapParse(arg, wordFreq);
  cout << "Remapping parse file took: " << difftime(time(NULL), start_wc)
       << " wall clock seconds\n";
  cout << "==== Elapsed time: " << difftime(time(NULL), start_main)
       << " wall clock seconds\n";
  return 0;
}
