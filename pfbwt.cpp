/* **************************************************************************
 * pfbwt.cpp
 * Output the BWT, the SA (option -S) or the sampled SA (option -s)
 * computed using the prefix-free parsing technique
 *
 * Usage:
 *   pfbwt[NT][64].x -h
 * for usage info
 *
 * See newscan.cpp for a description of what it does
 *
 **************************************************************************** */
#include <algorithm>
#include <assert.h>
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
#include "gsa/gsacak.h"
#include "utils.h"
}

using namespace std;
using namespace __gnu_cxx;

// -------------------------------------------------------------
// struct containing command line parameters and other globals
struct Args {
  char *basename;
  string parseExt = EXTPARSE; // extension final parse file
  string occExt = EXTOCC;     // extension occurrences file
  string dictExt = EXTDICT;   // extension dictionary file
  string lastExt = EXTLST;    // extension file containing last chars
  string saExt = EXTSAI;      // extension file containing sa info
  int w = 10;                 // sliding window size and its default
  int th = 0;                 // number of helper threads, default none
  bool SA = false;            // output all SA values
  int sampledSA = 0;          // output sampled SA values
};

// mask for sampled SA: start of a BWT run, end of a BWT run, or both
#define START_RUN 1
#define END_RUN 2

static long get_num_words(uint8_t *d, long n);
static long binsearch(uint_t x, uint_t a[], long n);
static int_t getlen(uint_t p, uint_t eos[], long n, uint32_t *seqid);
static void compute_dict_bwt_lcp(uint8_t *d, long dsize, long dwords, int w,
                                 uint_t **sap, int_t **lcpp);
static void fwrite_chars_same_suffix(vector<uint32_t> &id2merge,
                                     vector<uint8_t> &char2write,
                                     uint32_t *ilist, uint32_t *istart,
                                     FILE *fbwt, long &easy_bwts,
                                     long &hard_bwts);
static void fwrite_chars_same_suffix_sa(vector<uint32_t> &id2merge,
                                        vector<uint8_t> &char2write,
                                        uint32_t *ilist, uint32_t *istart,
                                        FILE *fbwt, long &easy_bwts,
                                        long &hard_bwts, int_t suffixLen,
                                        FILE *safile, uint8_t *bwsainfo, long);
static void fwrite_chars_same_suffix_ssa(
    vector<uint32_t> &id2merge, vector<uint8_t> &char2write, uint32_t *ilist,
    uint32_t *istart, FILE *fbwt, long &easy_bwts, long &hard_bwts,
    int_t suffixLen, FILE *ssafile, FILE *esafile, uint8_t *bwsainfo, long,
    int &, uint64_t &, int);
static uint8_t *load_bwsa_info(Args &arg, long n);

// class representing the suffix of a dictionary word
// instances of this class are stored to a heap to handle the hard bwts
struct SeqId {
  uint32_t id;   // lex. id of the dictionary word to which the suffix belongs
  int remaining; // remaining copies of the suffix to be considered
  uint32_t *bwtpos;   // list of bwt positions of this dictionary word
  uint8_t char2write; // char to be written (is the one preceeding the suffix)

  // constructor
  SeqId(uint32_t i, int r, uint32_t *b, int8_t c)
      : id(i), remaining(r), bwtpos(b) {
    char2write = c;
  }

  // advance to the next bwt position, return false if there are no more
  // positions
  bool next() {
    remaining--;
    bwtpos += 1;
    return remaining > 0;
  }
  bool operator<(const SeqId &a);
};

uint8_t last_char(Args &arg, uint32_t word, vector<string> dict_word) {
  uint8_t res = dict_word[word][dict_word[word].length() - arg.w - 1];
  return res;
}

bool SeqId::operator<(const SeqId &a) { return *bwtpos > *(a.bwtpos); }

/* *******************************************************************
 * Computation of the final BWT
 *
 * istart[] and islist[] are used together. For each dictionary word i
 * (considered in lexicographic order) for k=istart[i]...istart[i+1]-1
 * ilist[k] contains the ordered positions in BWT(P) containing word i
 * ******************************************************************* */
void bwt(Args &arg, uint8_t *d, long dsize, // dictionary and its size
         uint32_t *ilist, vector<vector<uint32_t>> children,
         long psize, // ilist, last and their size
         uint32_t *istart,
         long dwords, // starting point in ilist for each word and # words*
         vector<string> dict_word) // dictionary organized by work to make it
                                   // easier so get last
{
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
  // set d[0]==0 as this is the EOF char in the final BWT
  // assert(d[0] == Dollar);
  // d[0] = 0;

  // derive eos from sa. for i=0...dwords-1, eos[i] is the eos position of
  // string i in d
  uint_t *eos = sa + 1;
  for (int i = 0; i < dwords - 1; i++)
    assert(eos[i] < eos[i + 1]);

  // open output file
  FILE *fbwt = open_aux_file(arg.basename, "bwt", "wb");

  // main loop: consider each entry in the SA of dict
  time_t start = time(NULL);
  long full_words = 0;
  long easy_bwts = 0;
  long hard_bwts = 0;
  long next;
  uint32_t seqid;
  int lastbwt = Dollar;         // this is certainly not a BWT char
  uint64_t lastSa = UINT64_MAX; // this is an invalid SA entry

  // Adding the first characters for wich the contest is only w dollars
  for (uint32_t w : children[0]) {
    // compute next bwt char
    int nextbwt = last_char(arg, w - 1, dict_word);
    // in any case output BWT char
    if (fputc(nextbwt, fbwt) == EOF)
      die("BWT write error 0");
    lastbwt = nextbwt; // update lastbwt
    easy_bwts++;
  }

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
    if (sa[i] == 0 || d[sa[i] - 1] == EndOfWord) {
      full_words++;
      for (uint32_t w : children[seqid + 1]) {
        // compute next bwt char
        uint8_t nextbwt = last_char(arg, w - 1, dict_word);
        // in any case output BWT char
        if (fputc(nextbwt, fbwt) == EOF)
          die("BWT write error 0");
        lastbwt = nextbwt; // update lastbwt
        easy_bwts++;
      }
      continue; // proceed with next i
    }
    // ----- hard case: there can be a group of equal suffixes starting at i
    // save seqid and the corresponding char
    vector<uint32_t> id2merge(1, seqid);
    vector<uint8_t> char2write(1, d[sa[i] - 1]);
    /*cout << "seqid, char: " << seqid << ", " << d[sa[i] - 1] << endl;
    cout << "next: " << next << endl;
    cout << "lcp[next]: " << lcp[next] << endl;*/
    while (next < dsize && lcp[next] >= suffixLen) {
      assert(lcp[next] ==
             suffixLen); // the lcp cannot be greater than suffixLen
      // sa[next] cannot be a full word
      assert(sa[next] > 0 && d[sa[next] - 1] != EndOfWord);
      int_t nextsuffixLen = getlen(sa[next], eos, dwords, &seqid);
      // cout << "nextsuffixLen: " << nextsuffixLen << endl;
      assert(nextsuffixLen >= suffixLen);
      if (nextsuffixLen == suffixLen) {
        // cout << "in hard case at char: " << easy_bwts + hard_bwts << endl;
        id2merge.push_back(seqid);             // sequence to consider
        char2write.push_back(d[sa[next] - 1]); // corresponding char
        next++;
      } else
        break;
    }
    // cout << "next after: " << next << endl;
    // output to fbwt the bwt chars corresponding to the current dictionary
    // suffix, and, if requested, some SA values
    if (arg.SA)
      fwrite_chars_same_suffix_sa(id2merge, char2write, ilist, istart, fbwt,
                                  easy_bwts, hard_bwts, suffixLen, safile,
                                  bwsainfo, psize);
    else if (arg.sampledSA != 0)
      fwrite_chars_same_suffix_ssa(id2merge, char2write, ilist, istart, fbwt,
                                   easy_bwts, hard_bwts, suffixLen, ssafile,
                                   esafile, bwsainfo, psize, lastbwt, lastSa,
                                   arg.sampledSA);
    else
      fwrite_chars_same_suffix(id2merge, char2write, ilist, istart, fbwt,
                               easy_bwts, hard_bwts);
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

void print_help(char **argv, Args &args) {
  cout << "Usage: " << argv[0] << " <input filename> [options]" << endl;
  cout << "  Options: " << endl
       << "\t-w W\tsliding window size, def. " << args.w << endl
#ifndef NOTHREADS
       << "\t-t M\tnumber of helper threads, def. none " << endl
#endif
       << "\t-h  \tshow help and exit" << endl
       << "\t-s  \tcompute sampled suffix array" << endl
       << "\t-S  \tcompute full suffix array" << endl;
  exit(1);
}

void parseArgs(int argc, char **argv, Args &arg) {
  int c;
  extern char *optarg;
  extern int optind;

  puts("==== Command line:");
  for (int i = 0; i < argc; i++)
    printf(" %s", argv[i]);
  puts("");

  string sarg;
  while ((c = getopt(argc, argv, "t:w:sehS")) != -1) {
    switch (c) {
    case 's':
      arg.sampledSA |= START_RUN;
      break; // record SA position for start of runs
    case 'e':
      arg.sampledSA |= END_RUN;
      break; // record SA position for end of runs
    case 'S':
      arg.SA = true;
      break;
    case 'w':
      sarg.assign(optarg);
      arg.w = stoi(sarg);
      break;
    case 't':
      sarg.assign(optarg);
      arg.th = stoi(sarg);
      break;
    case 'h':
      print_help(argv, arg);
      exit(1);
    case '?':
      cout << "Unknown option. Use -h for help." << endl;
      exit(1);
    }
  }
  // the only input parameter is the file name
  arg.basename = NULL;
  if (argc == optind + 1)
    arg.basename = argv[optind];
  else {
    cout << "Invalid number of arguments" << endl;
    print_help(argv, arg);
  }
  // check algorithm parameters
  if (arg.SA && arg.sampledSA != 0) {
    cout << "You can either require the sampled SA or the full SA, not both";
    exit(1);
  }
  if (arg.w < 4) {
    cout << "Windows size must be at least 4\n";
    exit(1);
  }
#ifdef NOTHREADS
  if (arg.th != 0) {
    cout << "The NT version cannot use threads\n";
    exit(1);
  }
#else
  if (arg.th < 0) {
    cout << "Number of threads cannot be negative\n";
    exit(1);
  }
#endif
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
  parseArgs(argc, argv, arg);
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

  // compute and write the final bwt
  bwt(arg, d, dsize, ilist, children, psize, occ, dwords, dict_word);

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
// concatenation of dictionary words in lex order with EndOfWord (0x1) at the
// end of each word and d[size-1] = EndOfDict (0x0) at the very end. It is
// d[0]=Dollar (0x2) since the first word starts with $. There is another word
// somewhere ending with Dollar^wEndOfWord (it is the last word in the parsing,
// but its lex rank is unknown).
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

  *sap = sa;
  *lcpp = lcp;
}

// write to the bwt all the characters preceding a given suffix
// doing a merge operation if necessary
static void fwrite_chars_same_suffix(vector<uint32_t> &id2merge,
                                     vector<uint8_t> &char2write,
                                     uint32_t *ilist, uint32_t *istart,
                                     FILE *fbwt, long &easy_bwts,
                                     long &hard_bwts) {
  size_t numwords =
      id2merge.size(); // numwords dictionary words contain the same suffix
  bool samechar = true;
  for (size_t i = 1; (i < numwords) && samechar; i++)
    samechar = (char2write[i - 1] == char2write[i]);
  if (samechar) {
    for (size_t i = 0; i < numwords; i++) {
      uint32_t s = id2merge[i];
      for (long j = istart[s]; j < istart[s + 1]; j++) {
        if (fputc(char2write[0], fbwt) == EOF)
          die("BWT write error 1");
        easy_bwts++;
      }
    }
  } else {              // many words, many chars...
    vector<SeqId> heap; // create heap
    for (size_t i = 0; i < numwords; i++) {
      uint32_t s = id2merge[i];
      heap.push_back(SeqId(s, istart[s + 1] - istart[s], ilist + istart[s],
                           char2write[i]));
    }
    std::make_heap(heap.begin(), heap.end());
    while (heap.size() > 0) {
      // output char for the top of the heap
      SeqId s = heap.front();
      if (fputc(s.char2write, fbwt) == EOF)
        die("BWT write error 2");
      hard_bwts += 1;
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

// write to the bwt all the characters preceding a given suffix
// and the corresponding SA entries doing a merge operation
static void fwrite_chars_same_suffix_sa(
    vector<uint32_t> &id2merge, vector<uint8_t> &char2write, uint32_t *ilist,
    uint32_t *istart, FILE *fbwt, long &easy_bwts, long &hard_bwts,
    int_t suffixLen, FILE *safile, uint8_t *bwsainfo, long n) {
  size_t numwords =
      id2merge.size(); // numwords dictionary words contain the same suffix
  if (numwords == 1) {
    uint32_t s = id2merge[0];
    for (long j = istart[s]; j < istart[s + 1]; j++) {
      if (fputc(char2write[0], fbwt) == EOF)
        die("BWT write error 1");
      uint64_t sa = get_myint(bwsainfo, n, ilist[j]) - suffixLen;
      if (fwrite(&sa, SABYTES, 1, safile) != 1)
        die("SA write error 1");
    }
    easy_bwts += istart[s + 1] - istart[s];
  } else {              // many words, many chars...
    vector<SeqId> heap; // create heap
    for (size_t i = 0; i < numwords; i++) {
      uint32_t s = id2merge[i];
      heap.push_back(SeqId(s, istart[s + 1] - istart[s], ilist + istart[s],
                           char2write[i]));
    }
    std::make_heap(heap.begin(), heap.end());
    while (heap.size() > 0) {
      // output char for the top of the heap
      SeqId s = heap.front();
      if (fputc(s.char2write, fbwt) == EOF)
        die("BWT write error 2");
      uint64_t sa = get_myint(bwsainfo, n, *(s.bwtpos)) - suffixLen;
      if (fwrite(&sa, SABYTES, 1, safile) != 1)
        die("SA write error 2");
      hard_bwts += 1;
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

// write to the bwt all the characters preceding a given suffix
// and the corresponding sampled SA entries doing a merge operation
static void fwrite_chars_same_suffix_ssa(
    vector<uint32_t> &id2merge, vector<uint8_t> &char2write, uint32_t *ilist,
    uint32_t *istart, FILE *fbwt, long &easy_bwts, long &hard_bwts,
    int_t suffixLen, FILE *ssafile, FILE *esafile, uint8_t *bwsainfo, long n,
    int &bwtlast, uint64_t &salast, int ssa) {
  size_t numwords =
      id2merge.size(); // numwords dictionary words contain the same suffix
  if (numwords == 1) {
    // there is a single run, so a single potential SA value
    uint32_t s = id2merge[0];
    int bwtnext = char2write[0];
    if (bwtnext != bwtlast) {
      uint64_t pos = easy_bwts + hard_bwts;
      if (ssa & START_RUN) {
        uint64_t sa = get_myint(bwsainfo, n, ilist[istart[s]]) - suffixLen;
        if (fwrite(&pos, SABYTES, 1, ssafile) != 1)
          die("sampled SA write error 1a");
        if (fwrite(&sa, SABYTES, 1, ssafile) != 1)
          die("sampled SA write error 1b");
      }
      if (ssa & END_RUN) {
        pos--;
        if (fwrite(&pos, SABYTES, 1, esafile) != 1)
          die("sampled SA write error 1c");
        if (fwrite(&salast, SABYTES, 1, esafile) != 1)
          die("sampled SA write error 1d");
      }
      bwtlast = bwtnext;
    }
    for (long j = istart[s]; j < istart[s + 1]; j++) // write all BWT chars
      if (fputc(bwtnext, fbwt) == EOF)
        die("BWT write error 1");
    easy_bwts += istart[s + 1] - istart[s];
    if (ssa & END_RUN) // save the last sa value
      salast = get_myint(bwsainfo, n, ilist[istart[s + 1] - 1]) - suffixLen;
  } else {              // many words, many chars...
    vector<SeqId> heap; // create heap
    for (size_t i = 0; i < numwords; i++) {
      uint32_t s = id2merge[i];
      heap.push_back(SeqId(s, istart[s + 1] - istart[s], ilist + istart[s],
                           char2write[i]));
    }
    std::make_heap(heap.begin(), heap.end());
    while (heap.size() > 0) {
      // output char for the top of the heap
      uint64_t sa;
      SeqId s = heap.front();
      int bwtnext = s.char2write;
      if (fputc(bwtnext, fbwt) == EOF)
        die("BWT write error 2");
      if ((ssa & END_RUN) || (bwtnext != bwtlast))
        sa = get_myint(bwsainfo, n, *(s.bwtpos)) - suffixLen;
      if (bwtnext != bwtlast) {
        uint64_t pos = easy_bwts + hard_bwts;
        if (ssa & START_RUN) {
          if (fwrite(&pos, SABYTES, 1, ssafile) != 1)
            die("sampled SA write error 2a");
          if (fwrite(&sa, SABYTES, 1, ssafile) != 1)
            die("sampled SA write error 2b");
        }
        if (ssa & END_RUN) {
          pos--;
          if (fwrite(&pos, SABYTES, 1, esafile) != 1)
            die("sampled SA write error 2c");
          if (fwrite(&salast, SABYTES, 1, esafile) != 1)
            die("sampled SA write error 2d");
        }
        bwtlast = bwtnext;
      }
      if (ssa & END_RUN)
        salast = sa; // save current sa
      hard_bwts += 1;
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
