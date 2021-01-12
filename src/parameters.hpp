#include <iostream>
#include <string>

// file name extensions
#define EXTPARSE "parse"
#define EXTPARS0 "parse_old"
#define EXTOCC "occ"
#define EXTDICT "dict"
#define EXTDICZ "dicz"
#define EXTDZLEN "dicz.len"
#define EXTLST "last"
#define EXTCHILD "full_children"
#define EXTLIM "limits"
#define EXTBWTLIM "bwt_limits"
#define EXTSAI "sai"
#define EXTBWSAI "bwsai"
#define EXTILIST "ilist"
#define EXTSA "sa"
#define EXTSSA "ssa"
#define EXTESA "esa"

// mask for sampled SA: start of a BWT run, end of a BWT run, or both (for
// pfbwt_parseargs)
#define START_RUN 1
#define END_RUN 2

// -------------------------------------------------------------
// struct containing command line parameters and other globals for scan.cpp
struct Args {
  std::string inputFileName = "";
  std::string ReadsFileName = "";
  int w = 10;                      // sliding window size and its default
  int p = 100;                     // modulus for establishing stopping w-tuples
  bool SAinfo = false;             // compute SA information
  int th = 0;                      // number of helper threads
  int verbose = 0;                 // verbosity level
  std::string parseExt = EXTPARSE; // extension final parse file
  std::string occExt = EXTOCC;     // extension occurrences file
  std::string dictExt = EXTDICT;   // extension dictionary file
  std::string lastExt = EXTLST;    // extension file containing last chars
  std::string saExt = EXTSAI;      // extension file containing sa info
  char *basename;
  bool SA = false;   // output all SA values
  int sampledSA = 0; // output sampled SA values
  bool debug = false;
  bool rle = true;
};

void scan_parseargs(int argc, char **argv, Args &arg);
void pfbwt_parseargs(int argc, char **argv, Args &arg);
