#include "parameters.hpp"
#include <algorithm>
#include <assert.h>
#include <bits/stdint-uintn.h>
#include <ctime>
#include <errno.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <semaphore.h>
#include <sstream>
#include <stdexcept>
#include <stdint.h>
#include <string>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>
extern "C" {
#include "../external/gsa/gsacak.h"
#include "../external/utils.h"
}

using namespace std;
using namespace __gnu_cxx;

static long get_num_words(uint8_t *d, long n);
static long binsearch(uint_t x, uint_t a[], long n);
static int_t getlen(uint_t p, uint_t eos[], long n, uint32_t *seqid);
static void compute_dict_bwt_lcp(uint8_t *d, long dsize, long dwords, int w,
                                 uint_t **sap, int_t **lcpp);
static void fwrite_chars_same_suffix(
    vector<uint32_t> &id2merge, vector<uint8_t> &char2write, uint32_t *ilist,
    uint32_t *istart, FILE *fbwt, long &easy_bwts, long &hard_bwts,
    pair<uint32_t, uint32_t> *limits_ilist, vector<uint8_t> &pos);
static uint8_t *load_bwsa_info(Args &arg, long n);

// class representing the suffix of a dictionary word
// instances of this class are stored to a heap to handle the hard bwts
struct SeqId {
  uint32_t id;   // lex. id of the dictionary word to which the suffix belongs
  int remaining; // remaining copies of the suffix to be considered
  uint32_t *bwtpos; // list of bwt positions of this dictionary word
  pair<uint32_t, uint32_t> *limits_ilist; // list of the limits of this word
  uint8_t char2write; // char to be written (is the one preceeding the suffix)
  uint8_t pos;

  // constructor
  SeqId(uint32_t i, int r, uint32_t *b, int8_t c, pair<uint32_t, uint32_t> *d,
        uint8_t e)
      : id(i), remaining(r), bwtpos(b), limits_ilist(d), pos(e) {
    char2write = c;
  }

  // advance to the next bwt position, return false if there are no more
  // positions
  bool next() {
    remaining--;
    bwtpos += 1;
    limits_ilist++;
    return remaining > 0;
  }
  bool operator<(const SeqId &a);
  bool in_limits() {
    cout << "SedId in in_limits:" << limits_ilist->first
         << " <= " << unsigned(pos) << " < " << limits_ilist->second << " ?"
         << endl;
    return (pos < (int)limits_ilist->second && pos >= (int)limits_ilist->first);
  }
};

// Get the character at w+1 position from the end of the word.
uint8_t last_char(Args &arg, uint32_t word, vector<string> &dict_word) {
  uint8_t res = dict_word[word][dict_word[word].length() - arg.w - 1];
  return res;
}

bool SeqId::operator<(const SeqId &a) { return *bwtpos > *(a.bwtpos); }

/* *******************************************************************
 * Computation of the final BWT
 *
 * istart[] and islist[] are used together. For each dictionary word i
 * (considered in colexicographic order) for k=istart[i]...istart[i+1]-1
 * ilist[k] contains the ordered positions in BWT(P) containing word i
 * ******************************************************************* */
void bwt(Args &arg, uint8_t *d, long dsize, // dictionary and its size
         uint32_t *ilist, vector<vector<uint32_t>> children,
         long psize, // ilist, last and their size
         uint32_t *istart,
         long dwords, // starting point in ilist for each word and # words*
         vector<string> dict_word, // dictionary organized by word to make it
                                   // easier so get last
         pair<uint32_t, uint32_t> *limits_ilist) {
  // possibly read bwsa info file and open sa output file
  uint8_t *bwsainfo = load_bwsa_info(arg, psize);
  FILE *safile = NULL, *ssafile = NULL, *esafile = NULL;
  // open the necessary (sampled) SA files (possibly none)
  if (arg.SA)
    safile = open_aux_file(arg.basename, EXTSA, "wb");
  if (arg.sampledSA & START_RUN)
    ssafile = open_aux_file(arg.basename, EXTSSA, "wb");
  if (arg.sampledSA & END_RUN)
    esafile = open_aux_file(arg.basename, EXTESA, "wb");

  // compute sa and bwt of d and do some checking on them
  uint_t *sa;
  int_t *lcp;
  compute_dict_bwt_lcp(d, dsize, dwords, arg.w, &sa, &lcp);

  // derive eos from sa. for i=0...dwords-1, eos[i] is the eos position of
  // string i in d
  uint_t *eos = sa + 1;
  for (int i = 0; i < dwords - 1; i++)
    assert(eos[i] < eos[i + 1]);

  // open output file
  FILE *fbwt = open_aux_file(arg.basename, "bwt", "wb");

  cout << "Opening the output file" << endl;
  // main loop: consider each entry in the SA of dict
  time_t start = time(NULL);
  long full_words = 0;
  long easy_bwts = 0;
  long hard_bwts = 0;
  long next;
  uint32_t seqid;
  uint64_t lastSa = UINT64_MAX; // this is an invalid SA entry

  // Adding the first characters for wich the contest is only w dollars
  for (uint32_t w : children[0]) {
    cout << "w: " << w << endl;
    // compute next bwt char
    cout << "Compute nextbwt for last_char: " << w - 1 << endl;
    int nextbwt = last_char(arg, w - 1, dict_word);
    cout << "nextbwt: " << char(nextbwt) << endl;
    // in any case output BWT char
    if (fputc(nextbwt, fbwt) == EOF)
      die("BWT write error 0");
    easy_bwts++;
  }
  cout << "Finished adding the first characters" << endl;

  for (long i = dwords + 1; i < dsize; i = next) {
    // we are considering d[sa[i]....]
    next = i + 1; // prepare for next iteration
    // compute length of this suffix and sequence it belongs
    int_t suffixLen = getlen(sa[i], eos, dwords, &seqid);
    // ignore suffixes of lenght <= w
    if (suffixLen <= arg.w) {
      continue;
    }
    // ----- simple case: the suffix is a full word
    // ----- then we ouput the characters that folows in the tree
    if (sa[i] == 0 || d[sa[i] - 1] == EndOfWord) {
      full_words++;
      for (uint32_t w : children[seqid + 1]) {
        // compute next bwt char
        uint8_t nextbwt = last_char(arg, w - 1, dict_word);
        cout << "char written nextbwt: " << char(nextbwt) << endl;
        // in any case output BWT char
        if (fputc(nextbwt, fbwt) == EOF)
          die("BWT write error 0");
        easy_bwts++;
      }
      continue; // proceed with next i
    }
    // ----- hard case: there can be a group of equal suffixes starting at i
    // save seqid and the corresponding char
    vector<uint32_t> id2merge(1, seqid);
    vector<uint8_t> char2write(1, d[sa[i] - 1]);
    // vector<uint8_t> pos2test(1, dict_word[seqid].length() - suffixLen - 1);
    vector<uint8_t> pos2test(1, suffixLen);
    while (next < dsize && lcp[next] >= suffixLen &&
           (d[sa[next] - 1] != EndOfWord)) {
      cout << (next < dsize && lcp[next] >= suffixLen) << endl;
      // the lcp cannot be greater than suffixLen
      assert(lcp[next] == suffixLen);
      // sa[next] cannot be a full word
      // cout << sa[next] << endl;
      // cout << d[sa[next] - 1] << endl;
      // cout << (d[sa[next] - 1] == EndOfWord) << endl;
      assert(sa[next] > 0 && d[sa[next] - 1] != EndOfWord);
      int_t nextsuffixLen = getlen(sa[next], eos, dwords, &seqid);
      assert(nextsuffixLen >= suffixLen);
      if (nextsuffixLen == suffixLen) {
        id2merge.push_back(seqid);             // sequence to consider
        char2write.push_back(d[sa[next] - 1]); // corresponding char
        pos2test.push_back(suffixLen);
        // pos2test.push_back(dict_word[seqid].length() - suffixLen - 1);
        next++;
      } else
        break;
    }
    cout << "full_words: " << full_words << endl;
    cout << "easy_bwts: " << easy_bwts << endl;
    cout << "hard_bwts: " << hard_bwts << endl;
    // output to fbwt the bwt chars corresponding to the current dictionary
    // suffix
    fwrite_chars_same_suffix(id2merge, char2write, ilist, istart, fbwt,
                             easy_bwts, hard_bwts, limits_ilist, pos2test);
  }
  // write very last Sa pair
  if (arg.sampledSA & END_RUN) {
    uint64_t pos = easy_bwts + hard_bwts - 1;
    if (fwrite(&pos, SABYTES, 1, esafile) != 1)
      die("sampled SA write error 0e");
    if (fwrite(&lastSa, SABYTES, 1, esafile) != 1)
      die("sampled SA write error 0f");
  }
  assert(full_words == dwords);
  cout << "Full words: " << full_words << endl;
  cout << "Easy bwt chars: " << easy_bwts << endl;
  cout << "Hard bwt chars: " << hard_bwts << endl;
  cout << "Generating the final BWT took " << difftime(time(NULL), start)
       << " wall clock seconds\n";
  fclose(fbwt);
  delete[] lcp;
  delete[] sa;
  if (arg.SA or arg.sampledSA != 0)
    free(bwsainfo);
  if (arg.SA)
    fclose(safile);
  if (arg.sampledSA & 1)
    fclose(ssafile);
  if (arg.sampledSA & 2)
    fclose(esafile);
}

template <typename T = uint64_t> T read_binary64(ifstream &stream) {
  T a;
  stream.read((char *)&a, sizeof(a));
  return a;
}

int main(int argc, char **argv) {
  time_t start = time(NULL);

  // translate command line parameters
  Args arg;
  pfbwt_parseargs(argc, argv, arg);
  // read dictionary file
  FILE *g = open_aux_file(arg.basename, EXTDICT, "rb");
  fseek(g, 0, SEEK_END);
  long dsize = ftell(g);
  if (dsize < 0)
    die("ftell dictionary");
  if (dsize <= 1 + arg.w)
    die("invalid dictionary file");
  cout << "Dictionary file size: " << dsize << endl;
#if !M64
  if (dsize > 0x7FFFFFFE) {
    printf("Dictionary size greater than  2^31-2!\n");
    printf("Please use 64 bit version\n");
    exit(1);
  }
#endif

  uint8_t *d = new uint8_t[dsize];
  rewind(g);
  long e = fread(d, 1, dsize, g);
  if (e != dsize)
    die("fread");
  fclose(g);

  // Build a simpler dictionary without dollars
  vector<string> dict_word;
  string current_string = "";
  for (long i = 0; i < dsize; i++) {
    char c = d[i];
    if (c == EndOfWord) {
      dict_word.push_back(current_string);
      current_string = "";
    } else if (c == Dollar) {
      current_string += "$";
    } else if (c != EndOfDict) {
      current_string += c;
    }
  }

  // read occ file
  g = open_aux_file(arg.basename, EXTOCC, "rb");
  fseek(g, 0, SEEK_END);
  e = ftell(g);
  if (e < 0)
    die("ftell occ file");
  if (e % 4 != 0)
    die("invalid occ file");
  int dwords = e / 4;
  cout << "Dictionary words: " << dwords << endl;
  uint32_t *occ =
      new uint32_t[dwords + 1]; // dwords+1 since istart overwrites occ
  rewind(g);
  e = fread(occ, 4, dwords, g);
  if (e != dwords)
    die("fread 2");
  fclose(g);
  assert(dwords == get_num_words(d, dsize));

  // read ilist file
  g = open_aux_file(arg.basename, EXTILIST, "rb");
  fseek(g, 0, SEEK_END);
  e = ftell(g);
  if (e < 0)
    die("ftell ilist file");
  if (e % 4 != 0)
    die("invalid ilist file");
  long psize = e / 4;
  cout << "Parsing size: " << psize << endl;
  if (psize > 0xFFFFFFFEL)
    die("More than 2^32 -2 words in the parsing");
  uint32_t *ilist = new uint32_t[psize];
  rewind(g);
  e = fread(ilist, 4, psize, g);
  if (e != psize)
    die("fread 3");
  fclose(g);

  // convert occ entries into starting positions inside ilist
  // ilist also contains the position of EOF but we don't care about it since it
  // is not in dict
  // starting position in ilist of the smallest dictionary word
  uint32_t last = 1;
  for (long i = 0; i < dwords; i++) {
    uint32_t tmp = occ[i];
    occ[i] = last;
    last += tmp;
  }
  assert(last == psize);
  occ[dwords] = psize;

  // read children file
  ifstream file(string(arg.basename) + ".full_children", ios::in | ios::binary);
  uint32_t parsechar;
  file.read((char *)&parsechar, sizeof(parsechar));
  uint32_t sep = parsechar;
  vector<vector<uint32_t>> children(parsechar);
  uint32_t word = 0;
  while (file.read((char *)&parsechar, sizeof(parsechar))) {
    if (file.eof())
      break;
    if (parsechar == sep) {
      word++;
    } else {
      children[word].push_back(parsechar);
    }
  }

  // read the limits file
  vector<pair<uint32_t, uint32_t>> limits_ilist;
  ifstream ilist_limit_file(string(arg.basename) + ".limits_ilist",
                            ios::in | ios::binary);
  for (long i = 0; i < psize - 1; i++) {
    uint32_t l_start;
    ilist_limit_file.read((char *)&l_start, sizeof(l_start));
    uint32_t l_end;
    ilist_limit_file.read((char *)&l_end, sizeof(l_end));
    limits_ilist.push_back(make_pair(l_start, l_end));
    cout << "limits_ilist: " << l_start << " " << l_end << endl;
  }

  // compute and write the final bwt
  bwt(arg, d, dsize, ilist, children, psize, occ, dwords, dict_word,
      &limits_ilist[0]);

  delete[] ilist;
  delete[] occ;
  delete[] d;
  cout << "==== Elapsed time: " << difftime(time(NULL), start)
       << " wall clock seconds\n";
  return 0;
}

// --------------------- aux functions ----------------------------------

static uint8_t *load_bwsa_info(Args &arg, long n) {
  // maybe sa info is not really needed
  if (arg.SA == false and arg.sampledSA == 0)
    return NULL;
  // open .bwsa file for reading and .bwlast for writing
  FILE *fin = open_aux_file(arg.basename, EXTBWSAI, "rb");
  // allocate and load the bwsa array
  uint8_t *sai = (uint8_t *)malloc(n * IBYTES);
  if (sai == NULL)
    die("malloc failed (BWSA INFO)");
  long s = fread(sai, IBYTES, n, fin);
  if (s != n)
    die("bwsa info read");
  if (fclose(fin) != 0)
    die("bwsa info file close");
  return sai;
}

// compute the number of words in a dictionary
static long get_num_words(uint8_t *d, long n) {
  long i, num = 0;
  for (i = 0; i < n; i++)
    if (d[i] == EndOfWord)
      num++;
  assert(d[n - 1] == EndOfDict);
  return num;
}

// binary search for x in an array a[0..n-1] that doesn't contain x
// return the lowest position that is larger than x
static long binsearch(uint_t x, uint_t a[], long n) {
  long lo = 0;
  long hi = n - 1;
  while (hi > lo) {
    assert(((lo == 0) || x > a[lo - 1]) && x < a[hi]);
    int mid = (lo + hi) / 2;
    assert(x != a[mid]); // x is not in a[]
    if (x < a[mid])
      hi = mid;
    else
      lo = mid + 1;
  }
  assert(((hi == 0) || x > a[hi - 1]) && x < a[hi]);
  return hi;
}

// return the length of the suffix starting in position p.
// also write to seqid the id of the sequence containing that suffix
// n is the # of distinct words in the dictionary, hence the length of eos[]
static int_t getlen(uint_t p, uint_t eos[], long n, uint32_t *seqid) {
  assert(p < eos[n - 1]);
  *seqid = binsearch(p, eos, n);
  assert(eos[*seqid] > p); // distance between position p and the next $
  return eos[*seqid] - p;
}

// compute the SA and LCP array for the set of (unique) dictionary words
// using gSACA-K. Also do some checking based on the number and order of the
// special symbols d[0..dsize-1] is the dictionary consisting of the
// concatenation of reversed dictionary words in lex order (so the words of the
// dictionary are in colexicographic order. with EndOfWord (0x1) at the end of
// each word and d[size-1] = EndOfDict (0x0) at the very end. It is d[0]=Dollar
// (0x2) since the first word starts with $. There is another word somewhere
// ending with Dollar^wEndOfWord (it is the last word in the parsing, but its
// lex rank is unknown).
static void compute_dict_bwt_lcp(uint8_t *d, long dsize, long dwords, int w,
                                 uint_t **sap,
                                 int_t **lcpp) // output parameters
{
  uint_t *sa = new uint_t[dsize];
  int_t *lcp = new int_t[dsize];
  (void)dwords;
  (void)w;

  cout << "Each SA entry: " << sizeof(*sa) << " bytes\n";
  cout << "Each LCP entry: " << sizeof(*lcp) << " bytes\n";

  cout << "Computing SA and LCP of dictionary" << endl;
  time_t start = time(NULL);
  gsacak(d, sa, lcp, NULL, dsize);
  cout << "Computing SA/LCP took " << difftime(time(NULL), start)
       << " wall clock seconds\n";
  // ------ do some checking on the sa
  assert(d[dsize - 1] == EndOfDict);
  assert(sa[0] == (unsigned long)dsize - 1); // sa[0] is the EndOfDict symbol
  for (long i = 0; i < dwords; i++)
    assert(d[sa[i + 1]] == EndOfWord); // there are dwords EndOfWord symbols
  // EndOfWord symbols are in position order, so the last is d[dsize-2]
  assert(sa[dwords] == (unsigned long)dsize - 2);
  cout << "finished checking the SA" << endl;
  *sap = sa;
  *lcpp = lcp;
}

bool pos_in_limits(int_t pos, pair<uint32_t, uint32_t> word_limit) {
  cout << "pos_in_limits: " << word_limit.first << " <= " << unsigned(pos)
       << " < " << word_limit.second << " ?" << endl;
  return (pos < (int)word_limit.second && pos >= (int)word_limit.first);
}

// write to the bwt all the characters preceding a given suffix
// doing a merge operation if necessary
static void fwrite_chars_same_suffix(
    vector<uint32_t> &id2merge, vector<uint8_t> &char2write, uint32_t *ilist,
    uint32_t *istart, FILE *fbwt, long &easy_bwts, long &hard_bwts,
    pair<uint32_t, uint32_t> *limits_ilist, vector<uint8_t> &pos2test) {
  size_t numwords =
      id2merge.size(); // numwords dictionary words contain the same suffix
  bool samechar = true;
  for (size_t i = 1; (i < numwords) && samechar; i++)
    samechar = (char2write[i - 1] == char2write[i]);
  if (samechar) {
    for (size_t i = 0; i < numwords; i++) {
      uint32_t s = id2merge[i];
      // CHECK things around here
      auto limit_j = limits_ilist + istart[s] - 1;
      for (long j = istart[s]; j < istart[s + 1]; j++) {
        cout << "char to be written: " << char2write[0] << endl;
        if (pos_in_limits(pos2test[i], *limit_j)) {
          if (fputc(char2write[0], fbwt) == EOF)
            die("BWT write error 1");
          easy_bwts++;
        }
        limit_j++;
      }
    }
  } else {              // many words, many chars...
    vector<SeqId> heap; // create heap
    for (size_t i = 0; i < numwords; i++) {
      uint32_t s = id2merge[i];
      heap.push_back(SeqId(s, istart[s + 1] - istart[s], ilist + istart[s],
                           char2write[i], limits_ilist + istart[s],
                           pos2test[i]));
    }
    std::make_heap(heap.begin(), heap.end());
    while (heap.size() > 0) {
      // output char for the top of the heap
      SeqId s = heap.front();
      if (s.in_limits()) {
        if (fputc(s.char2write, fbwt) == EOF)
          die("BWT write error 2");
        hard_bwts += 1;
      }
      // remove top
      pop_heap(heap.begin(), heap.end());
      heap.pop_back();
      // if remaining positions, reinsert to heap
      if (s.next()) {
        heap.push_back(s);
        push_heap(heap.begin(), heap.end());
      }
    }
  }
}
