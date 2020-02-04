#include "prettyprint.hpp"
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

template <typename T = uint64_t> void write_binary(T item, ofstream &stream) {
  stream.write((char *)&item, sizeof(item));
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
graph_structure_t graph_structure(const string &filename,
                                  uint32_t &alphabet_parse) {
  graph_structure_t structure;
  ifstream file(filename + ".parse", ios::in | ios::binary);
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
  while (file.read((char *)&parsechar, sizeof(parsechar)),
         parsechar != separator) {
    // reading genome
    // cout << "character_in_the_parse " << parsechar << endl;
    if (not file) {
      cerr << "Bad file structure: " << filename + ".parse" << endl;
      return structure;
    }
    structure.push_back({parsechar, i, false});
    i++;
  }

  // genome is zero based 0 1 2 3 4 ... i
  uint64_t genome_length = i;
  get<2>(structure.back()) = true; // the last one is a leaf
  cout << "genome_parse_length: " << genome_length << endl;
  i++; // final node of the genome

  // Parsing the reads
  for (;;) {
    if (file.eof())
      break;
    uint64_t pos_in_genome = read_binary64(file);
    assert(pos_in_genome <= genome_length);
    uint32_t parsechar = read_binary(file);
    structure.push_back(
        {parsechar, pos_in_genome, false}); // link with the genome
    while (parsechar = read_binary(file), parsechar != separator) {
      if (file.eof())
        break;
      // cout << parsechar << endl;
      structure.push_back({parsechar, i, false});
      i++;
    }
    get<2>(structure.back()) = true; // the last one is a leaf
    // cout << "pos_in_genome, genome length: " << pos_in_genome << " "
    //     << genome_length << endl;
    i++;
    if (file.eof())
      break;
  }

  return structure;
}

using triple = tuple<uint64_t, uint64_t, uint64_t>;
pair<uint64_t, uint64_t> getPairofParsechar(triple &l, vector<triple> &in) {
  return pair<uint64_t, uint64_t>{get<1>(l), get<1>(in[get<0>(l)])};
}

#define all(x) begin(x), end(x)
vector<uint64_t> doubling_algorithm(graph_structure_t &graph_structure) {
  uint64_t n = graph_structure.size();
  vector<uint64_t> result(n + 1);
  vector<triple> intermediate;
  // At the begining this contains the triple (parent, parse_char, self)
  // At the end it should contain (0, SA index, self)
  vector<uint64_t> order(n + 1, 0);
  intermediate.push_back({0, 0, 0}); // Add the root
  uint64_t i = 1;
  for (auto &e : graph_structure) {
    // Fill the array with the indexes of the nodes to organize !
    intermediate.push_back({get<1>(e), get<0>(e), i});
    order[i] = i;
    i++;
  }
  // cout << intermediate << endl;
  // cout << order << endl;

  while (true) {
    stable_sort(all(order), [&](const uint64_t lhs, const uint64_t rhs) {
      auto &l = intermediate[lhs];
      auto &r = intermediate[rhs];
      return getPairofParsechar(l, intermediate) <
             getPairofParsechar(r, intermediate);
    });
    pair<uint64_t, uint64_t> prevpair = {0, 0};
    uint64_t crank = 0;
    bool cont = 0;
    for (uint64_t i = 0; i < n + 1; i++) {
      auto &el = intermediate[order[i]];
      auto p = getPairofParsechar(el, intermediate);
      if (prevpair != p) {
        crank++;
        prevpair = p;
      }
      get<1>(el) = crank;
      get<0>(el) = get<0>(intermediate[get<0>(el)]);
      if (get<0>(el) > 0) {
        cont = 1;
      }
    }
    if (!cont) {
      break;
    }
    // cout << intermediate << endl;
  }
  // To avoid equalities, put a total ordering
  i = 0;
  for (auto &c : order) {
    get<1>(intermediate[c]) = i;
    i++;
  }
  // Creating the Inverted suffix array
  for (auto &e : intermediate) {
    result[get<1>(e)] = get<2>(e); // 2 is self and 1 ordering in the SA
    assert(get<0>(e) == 0); // prefix_sort was succesful, for all node we went
                            // back until the root
  }
  return result;
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
  // read parse file
  uint32_t alphabet_parse = 0;
  auto structure = graph_structure(filename, alphabet_parse);
  // cout << "structure: " << structure << endl;
  cout << "Compute wheeler order and then SA via doubling algorithm" << endl;
  vector<uint64_t> sa = doubling_algorithm(structure);
  // cout << "sa: " << sa << endl;
  // inversing the graph so we can build a BWT
  cout << "Building the reverse graph for the BWT" << endl;
  uint64_t n = structure.size();
  vector<vector<uint32_t>> children_char(n + 1);
  for (auto &e : structure) {
    // fill with the chars of the children
    children_char[get<1>(e)].push_back(get<0>(e));
  }
  // cout << children_char << endl;
  vector<vector<uint32_t>> children_full_word(alphabet_parse + 1);
  children_full_word[0] = children_char[0];
  for (uint64_t i = 1; i < n; i++) {
    // link to the previous word
    if (get<1>(structure[i]) > 0) {
      auto e = structure[get<1>(structure[i]) - 1];
      children_full_word[get<0>(e)].insert(children_full_word[get<0>(e)].end(),
                                           children_char[i].begin(),
                                           children_char[i].end());
    }
  }
  // cout << children_full_word << endl;
  cout << "Compute the BWT from the SA" << endl;
  vector<uint32_t> BWT;
  BWT.push_back(0); // empty word
  for (uint64_t i = 0; i < n + 1; i++) {
    for (auto &c : children_char[sa[i]]) {
      BWT.push_back(c);
    }
  }
  // cout << "BWT: " << BWT << endl;
  // Creating the F vector
  ifstream file_occ(filename + ".occ", ios::in | ios::binary);
  vector<uint32_t> occ(alphabet_parse + 1, 0);
  for (uint32_t i = 1; i < alphabet_parse + 1; i++)
    occ[i] = read_binary(file_occ);
  occ[0] = 1; // empty word
  // cout << "occ: " << occ << endl;
  vector<uint32_t> F(alphabet_parse + 1, 0);
  for (uint32_t i = 1; i <= alphabet_parse; i++)
    F[i] = F[i - 1] + occ[i - 1];
  // cout << "F: " << F << endl;
  vector<uint32_t> ilist(n + 1, 0);
  for (uint32_t i = 0; i <= n; i++) {
    ilist[F[BWT[i]]++] = i;
    occ[BWT[i]]--;
  }
  // cout << "ilist: " << ilist << endl;

  for (uint32_t i = 0; i < alphabet_parse + 1; i++) {
    assert(occ[i] == 0);
  }

  cout << "Saving .full_children and .ilist files" << endl;
  // Saving children ordered by sa
  auto children_file = ofstream(filename + ".full_children");
  uint32_t children_sep = alphabet_parse + 1;
  for (uint32_t i = 0; i < alphabet_parse + 1; i++) {
    write_binary(children_sep, children_file);
    for (auto &c : children_full_word[i]) {
      write_binary(c, children_file);
    }
  }
  children_file.close();
  cout << ".full_children file writen and closed" << endl;
  // Saving ilist
  auto ilist_file = ofstream(filename + ".ilist");
  for (uint32_t i = 0; i < n + 1; i++)
    write_binary(ilist[i], ilist_file);
  ilist_file.close();
  cout << ".ilist file writen and closed" << endl;
  printf("==== Elapsed time: %.0lf wall clock seconds\n",
         difftime(time(NULL), start_wc));
  return 0;
}
