#include "parameters.hpp"
#include <unistd.h>

using namespace std;

void scan_print_help(char **argv, Args &args) {
  cout << "Usage: " << argv[0] << " <input filename> [options]" << endl;
  cout << "  Options: " << endl
       << "\t-w W\tsliding window size, def. " << args.w << endl
       << "\t-p M\tmodulo for defining phrases, def. " << args.p << endl
       << "\t-c  \tdiscard redundant information" << endl
       << "\t-h  \tshow help and exit" << endl
       << "\t-s  \tcompute suffix array info" << endl;
  exit(1);
}

void scan_parseargs(int argc, char **argv, Args &arg) {
  int c;
  extern char *optarg;
  extern int optind;

  puts("==== Command line:");
  for (int i = 0; i < argc; i++)
    printf(" %s", argv[i]);
  puts("");

  string sarg;
  while ((c = getopt(argc, argv, "p:w:sht:vc")) != -1) {
    switch (c) {
    case 's':
      arg.SAinfo = true;
      break;
    case 'w':
      sarg.assign(optarg);
      arg.w = stoi(sarg);
      break;
    case 'p':
      sarg.assign(optarg);
      arg.p = stoi(sarg);
      break;
    case 't':
      sarg.assign(optarg);
      arg.th = stoi(sarg);
      break;
    case 'v':
      arg.verbose++;
      break;
    case 'h':
      scan_print_help(argv, arg);
      exit(1);
    case '?':
      cout << "Unknown option. Use -h for help." << endl;
      exit(1);
    }
  }
  // the only input parameter is the file name
  if (argc == optind + 2) {
    arg.inputFileName.assign(argv[optind]);
    arg.ReadsFileName.assign(argv[optind + 1]);
  } else {
    cout << "Invalid number of arguments" << endl;
    scan_print_help(argv, arg);
  }
  // check algorithm parameters
  if (arg.w < 4) {
    cout << "Windows size must be at least 4\n";
    exit(1);
  }
  if (arg.p < 10) {
    cout << "Modulus must be at leas 10\n";
    exit(1);
  }
  if (arg.th != 0) {
    cout << "The NT version cannot use threads\n";
    exit(1);
  }
}

void pfbwt_print_help(char **argv, Args &args) {
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

void pfbwt_parseargs(int argc, char **argv, Args &arg) {
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
      pfbwt_print_help(argv, arg);
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
    pfbwt_print_help(argv, arg);
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
