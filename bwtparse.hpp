//
// Created by Jan Studený on 15/01/2020.
//
/*
  - Input
    - file.parse
    - file.dic
  - Intermediate
    - "file.diclast"
  - Output
    - file.ilist
    - file.occ
    - file.bwlast
- Steps:
  - 0. create .diclast from .dic
    - diclast[i] = dic[i][-w]
  - 1. from file.parse create a set of triples
    - the ids can be 0-k_0-1 (genome) and then one by one reads (k_0...k_1-1)
(...) (no need to be sorted anyhow because we need to be consistent within the
read only)
    - [(metachar,parent,is_leaf)]
  - 2. doubling algorithm on a set of "triples"
    - [(metachar,parent,is_leaf)] -> inverted "eXtended" suffix array (a.k.a.
wheeler order)
  - 3. for each word just output the position using the result from 2. (+ in
parallell bwlast as for bwlast[parent] = diclast[me] & "dollars/special
character" for uninitialised, end of reads)
    - vector<vector<int>> multiplicities (for every word where it occurs)
    - for each word, sort it
  - 3.5 put to old format
    -  multiplicites to .ilist .occ
    - .occ[i] = multiplicities[i].size()
  - 4. (external) run the pfbwt using the fake .ilist .occ .bwlast
 */

#include "prettyprint.hpp"
#include "utils.h"
#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#define DebugDollar 0x0
using namespace std;

#define all(x) begin(x), end(x)
template <typename T = uint32_t> using carr = pair<T *, long>;

using triple = tuple<int64_t, int64_t, int64_t>;
using pii = pair<int64_t, int64_t>;

pair<int64_t, int64_t> getPair(tuple<int64_t, int64_t, int64_t> &l,
                               vector<triple> &in) {
  return pair<int64_t, int64_t>{get<1>(l), get<1>(in[get<0>(l)])};
}


template <typename T = int64_t> T read_binary(ifstream &stream) {
  T a;
  stream.read((char *)&a, sizeof(a));
  return a;
}

template <> int64_t read_binary<int64_t>(ifstream &stream) {
  int64_t a;
  stream.read((char *)&a, sizeof(a));
  return __builtin_bswap64(a);
}

template <> int32_t read_binary<int32_t>(ifstream &stream) {
  int32_t a;
  stream.read((char *)&a, sizeof(a));
  return __builtin_bswap32(a);
}

template <typename T = int64_t> void write_binary(T item, ofstream &stream) {
  stream.write((char *)&item, sizeof(item));
}


string get_content(const string &filename) {
  std::ifstream t(filename);
  std::stringstream buffer;
  buffer << t.rdbuf();
  return buffer.str();
}

void save_content(const string &content, string filename) {
  std::ofstream t(filename);
  t << content;
  t.close();
}

template <typename T> void save_carr(T carr, ofstream &stream) {
  for (int64_t i = 0; i < carr.second; i++) {
    write_binary(carr.first[i], stream);
  }
}
