#include <assert.h>
#include <ctime>
#include <errno.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <sstream>
#include <stdexcept>
#include <stdint.h>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>
#include "prettyprint.hpp"
extern "C" {
#include "utils.h"
}

using namespace std;
using graph_structure_t = vector<tuple<uint32_t, uint64_t, bool>>;

template <typename T = uint32_t> T read_binary(ifstream &stream) {
  T a;
  stream.read((char *)&a, sizeof(a));
  return a;
}

template <typename T = uint64_t> T read_binary64(ifstream &stream) {
  T a;
  stream.read((char *)&a, sizeof(a));
  return a;
}

// structure: (character_in_the_parse, parent_in_the_tree, is_leaf)
/*
 * A first ordering of the structure
 *
 *          i i+1 i+2 i+3
 *         /
 *  0 1 2 3 4 5 6 7 8 ... i-1
 *
 *
 */
graph_structure_t graph_structure(const string &filename) {
    graph_structure_t structure;
    ifstream file(filename + ".parse", ios::in | ios::binary);
    uint32_t alphabet_parse = 0;
    file.read((char *)&alphabet_parse, sizeof(alphabet_parse));
    cout << "Size of the parse's alphabet: " << alphabet_parse << endl;
    uint32_t separator = alphabet_parse + 1;
    if (not file) {
      cerr << "Bad file structure: " << filename + ".parse" << endl;
      return structure;
    }

    // Parsing the genome
    int i = 0; // the root doesn't have a parent
    uint32_t parsechar;
    while (file.read((char *)&parsechar, sizeof(parsechar)), parsechar != separator) {
      // reading genome
      cout << "character_in_the_parse " << parsechar << endl;
      if (not file) {
        cerr << "Bad file structure: " << filename + ".parse" << endl;
        return structure;
      }
      structure.push_back({parsechar, i, false});
      i++;
    }

    // genome is zero based 0 1 2 3 4 ... i-1
    uint64_t genome_length = i;
    get<2>(structure.back()) = true; // the last one is a leaf
    cout << "genome_length: " << genome_length << endl;

    //Parsing the reads
    for (;;) {
      uint64_t pos_in_genome = read_binary64(file);
      uint32_t parsechar = read_binary(file);
      structure.push_back({parsechar, pos_in_genome, false}); //link with the genome
      while (parsechar = read_binary(file), parsechar != separator) {
        if (parsechar == 0) break; // alphabet is starting at 1
        structure.push_back({parsechar, i, false});
        pos_in_genome++;
        i++;
      }
      get<2>(structure.back()) = true; // the last one is a leaf
      assert(pos_in_genome <= genome_length);
      if (parsechar == 0) break; // alphabet is starting at 1
    }

    return structure;
}


int main(int argc, char *argv[]) {
  assert(argc == 2);
  string filename = argv[1];
  puts("==== Command line:");
  for (int i = 0; i < argc; i++)
    printf(" %s", argv[i]);
  puts("");
  // start measuring wall clock time
  time_t start_wc = time(NULL);
  cout << "Get the structure of the tree" << endl;
  auto structure = graph_structure(filename);
  cout << "structure: " << structure << endl;
  // read parse file
  printf("==== Elapsed time: %.0lf wall clock seconds\n", difftime(time(NULL),start_wc));
  return 0;
}
