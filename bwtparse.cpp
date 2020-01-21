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
    - the ids can be 0-k_0-1 (genome) and then one by one reads (k_0...k_1-1) (...) (no need to be sorted anyhow because we need to be consistent within the read only)
    - [(metachar,parent,is_leaf)]
  - 2. doubling algorithm on a set of "triples"
    - [(metachar,parent,is_leaf)] -> inverted "eXtended" suffix array (a.k.a. wheeler order)
  - 3. for each word just output the position using the result from 2. (+ in parallell bwlast as for bwlast[parent] = diclast[me] & "dollars/special character" for uninitialised, end of reads)
    - vector<vector<int>> multiplicities (for every word where it occurs)
    - for each word, sort it
  - 3.5 put to old format
    -  multiplicites to .ilist .occ
    - .occ[i] = multiplicities[i].size()
  - 4. (external) run the pfbwt using the fake .ilist .occ .bwlast
 */

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include "utils.h"
#include "prettyprint.hpp"

#define DebugDollar 0x0
using namespace std;


template <typename T=int64_t>
T read_binary(ifstream& stream)  {
    T a;
    stream.read((char*)&a, sizeof(a));
    return a;
}

template <>
int64_t read_binary<int64_t>(ifstream& stream)  {
    int64_t a;
    stream.read((char*)&a, sizeof(a));
    return __builtin_bswap64(a);
}

template <>
int32_t read_binary<int32_t>(ifstream& stream)  {
    int32_t a;
    stream.read((char*)&a, sizeof(a));
    return __builtin_bswap32(a);
}

template <typename T=int64_t>
void write_binary(T item,ofstream& stream)  {
    stream.write((char*)&item, sizeof(item));
}

/* [(metachar,parent,is_leaf)] */
/*
 *   i+1 i+2 i+3
 *              \
 * 1 2 3 4 5 6 7 8 ... i 0
 *
 *
 */
using graph_structure_t = vector<tuple<int64_t,int64_t,bool>>;
graph_structure_t graph_structure(const string& filename) {
    graph_structure_t structure;
    auto parse_filename = filename + ".parse";
    ifstream file (parse_filename, ios::in | ios::binary);
    //vector<int64_t> genome;
    const int64_t root = 0;
    structure.push_back({DebugDollar,root,0}); // dummy dollar in the back
    int i=1;
    int genome_length = 0;
    if (not file) {
        cerr << "Bad file" << parse_filename << endl;
        return structure;
    }
    int64_t metachar;
    bool is_leaf = false; // genome has "no leaves" (we want exactly one dollar)
    while(metachar = read_binary(file),metachar != PRIME+1) {
        // reading genome
        if (not file) {
            cerr << "Bad file structure: " << parse_filename << endl;
            return structure;
        }
        structure.push_back({metachar,i+1,is_leaf});
        is_leaf = false;
        i++;
    }
    genome_length = i;
    get<1>(structure.back()) = root; // dummy as last
    for(;;) {
        uint64_t pos_in_genome = read_binary(file); // I should take its last position
        if (not file) {break;}
        int64_t metachar;
        bool is_leaf = true;
        while(metachar = read_binary(file),metachar != PRIME+1) {
            structure.push_back({metachar,i+1,is_leaf});
            pos_in_genome++; // compensate that I have initial position instead of last position: FIXME
            i++;
            is_leaf = false;
        }
        // genome is one based now 1 2 3 4 ... i 0
        //pos_in_genome--; // one past end -> end
        assert(pos_in_genome < genome_length);
        get<1>(structure.back()) = pos_in_genome+1 == genome_length ? root : pos_in_genome + 1; // joining the genome
        // TODO: when everything is reversed make sure the pos_in_genome is also reversed (and correct)
    }
    return structure;
}
// TODO: change to structs
#define metachar 0
#define parent 1
#define place_in_xbwt second
#define all(x) begin(x),end(x)
template <typename T=uint32_t>
using carr = pair<T*,long>;

using triple = tuple<int64_t,int64_t,int64_t>;
using pii = pair<int64_t,int64_t>;

pair<int64_t, int64_t> getPair(tuple<int64_t, int64_t, int64_t> &l, vector<triple> &in) { return pair<int64_t,int64_t>{get<1>(l), get<1>(in[get<0>(l)])}; }

template <typename T>
void save_carr(T carr, ofstream& stream) {
    for(int64_t i=0; i < carr.second;i++) {
        write_binary(carr.first[i],stream);
    }
}

void prefix_sort(vector<triple>& in_out, bool total_order, graph_structure_t& graph_structure) {
    // [parent,metachar,me] -> [root,inverted_suffix_arry_idx,me]
    // assumptions (get<2> is a permutation from 0 to n-1)
    // assumptions (get<0> < get<2>)
    // time complexity O(nlog^2(n)) which can be easily O(nlogn)
    int64_t n = in_out.size();
    // sort them according to their
    stable_sort(all(in_out), [](const triple& lhs, const triple& rhs) {
        return get<2>(lhs) < get<2>(rhs);
    });
    vector<int64_t> order(n,0);
    for(int64_t i=0;i<n;i++) order[i] = i;
    while(true) {
        stable_sort(all(order),[&](const int64_t lhs,const int64_t rhs){
            auto& l = in_out[lhs];
            auto& r = in_out[rhs];
            return getPair(l, in_out) < getPair(r, in_out);
        });
        pair<int64_t,int64_t> prevpair = {-1,-1};
        int64_t crank = -1;
        bool cont = 0;
        for(int64_t i=0;i<n;i++) {
            auto& el = in_out[order[i]];
            auto p = getPair(el, in_out);
            if (prevpair != p) {
                crank++;
                prevpair = p;
            }
            get<1>(el) = crank;
            get<0>(el) = get<0>(in_out[get<0>(el)]);
            if (get<0>(el) > 0) cont = 1;
        }
        if (!cont) {
            break;
        }
        cout << in_out << endl;
    }
    if (total_order) {
        stable_sort(all(order),[&](const int64_t lhs,const int64_t rhs) // should be unnecessary
        {
            auto& l = in_out[lhs];
            auto& r = in_out[rhs];
            return pii{get<1>(l),get<2>(l)} < pii{get<1>(r),get<2>(r)};
        });
        int64_t i=0;
        for(auto& c : order) {
            auto is_leaf = get<2>(graph_structure[c]);
            get<1>(in_out[c]) = is_leaf ? -1 : i;
            i = is_leaf ? i : i + 1;
        }
    }
}

/* [(metachar,parent,is_leaf)] -> inverted "eXtended" suffix array (a.k.a. wheeler order) */
vector<int64_t> doubling_algorithm(graph_structure_t& graph_structure) {
    vector<triple> intermediate;
    vector<int64_t> result(graph_structure.size());
    int i=0;
    for(auto& e : graph_structure) {
        intermediate.push_back({get<parent>(e),get<metachar>(e),i});
        i++;
    }
    cout << intermediate << endl;
    prefix_sort(intermediate,true,graph_structure);
    for(auto& e : intermediate) {
        result[get<2>(e)] = get<1>(e);
        assert(get<0>(e) == 0); // prefix_sort was succesful
    }
    return result;
}


pair<vector<vector<int64_t>>, string>
fake_lists(vector<int64_t> &doubling_result, graph_structure_t &graph_structure,
           string &diclast, int64_t dictionary_size) {
    vector<vector<int64_t>> multiplicities(dictionary_size+1); // + 1 for dollar word (smallest word)
    string bwlast = string(diclast.size(), Dollar); // TODO: Check if dollar works
    for(int me=0;me<graph_structure.size();me++) {
        auto& node = graph_structure[me];
        multiplicities[get<metachar>(node)].push_back(doubling_result[get<parent>(node)]);
        bwlast[get<parent>(node)] = diclast[get<metachar>(node)]; // check for correctness
    }
    multiplicities[0][0] = doubling_result[1]; // the only dollar placed is from the genome
    for(auto& e : multiplicities) {
        sort(all(e));
    }
    return {multiplicities,bwlast};
}

// multiplicites to occurences
carr<> mults_to_occ(vector<vector<int64_t>> multiplicities) {
    // excluding the first (dollar) symbol
    auto dwords = multiplicities.size()-1;
    auto occ = new uint32_t[dwords];
    for(int i=1;i<multiplicities.size();i++) {
        occ[i-1] = multiplicities[i].size();
    }
    return {occ,dwords};
}

// multiplicity list transformation
pair<carr<>,carr<>> mult_list_to_ilist_istart(const vector<vector<int64_t>>& mult) {
    int dwords = mult.size();
    int elems = 0;
    for(auto& e : mult) {
        elems += e.size();
    }
    auto ilist = new uint32_t[elems];
    auto istart = new uint32_t[dwords];

    int i = 0;
    int j = 0;
    for(auto& e : mult) {
        istart[i] = j;
        for(auto& c : e) {
            ilist[j] = c;
            j++;
        }
        i++;
    }
    return {{ilist,elems},{istart,dwords}};
}

carr<uint8_t> dict_to_old_repr(const vector<string>& dict) {
    int tot_length = 0;
    for(auto& w : dict) {
        tot_length += w.size();
        tot_length++;
    }
    tot_length++;
    auto* res = new uint8_t[tot_length];
    int idx = 0;
    for(auto& w : dict) {
        for(auto& c : w) {
            res[idx] = c;
            idx++;
        }
        res[idx] = 0x01;
        idx++;
    }
    res[idx] = 0;
    return {res,tot_length};
}

vector<string> dict_from_old_repr(const string& repr) {
    vector<string> dict;
    int i=0;
    stringstream temp_str;
    for(auto& c : repr) {
        if (c == 0x1) {
            dict.push_back(temp_str.str());
            temp_str.str("");
        }
        else if (c == 0x0) {
            assert(i+1==repr.length());
            assert(temp_str.str()==""s);
        }
        else {
            temp_str << c;
        }
        i++;
    }
    return dict;
}

string get_content(const string& filename) {
    std::ifstream t(filename);
    std::stringstream buffer;
    buffer << t.rdbuf();
    return buffer.str();
}

void save_content(const string& content,string filename) {
    std::ofstream t(filename);
    t << content;
    t.close();
}

string get_diclast(const string& filename,int64_t window_size) {
    auto content = get_content(filename + ".dict");
    auto dictionary = dict_from_old_repr(content);
    cout << dictionary << endl;
    string diclast;
    for(auto& word : dictionary) {
        assert(word.length() > window_size);
        diclast.push_back(word[word.length()-window_size-1]); // test -1
    }
    return diclast;
}



int main(int argc, char** argv) {
    assert(argc==2);
    string filename = argv[1];
    int64_t window_size = 0;
    auto structure = graph_structure(filename);
    cout << structure << endl;
    auto doubling_result = doubling_algorithm(structure);
    cout << doubling_result << endl;
    auto diclast = get_diclast(filename,window_size);
    auto [multiplicities, bwlast] = fake_lists(doubling_result,structure,diclast,diclast.size());
    cout << multiplicities << endl;
    auto occ = mults_to_occ(multiplicities);
    cout << occ << endl;
    auto [ilist, istart] = mult_list_to_ilist_istart(multiplicities);
    cout << ilist << endl;
    cout << istart << endl;
    auto ilist_file = ofstream(filename + ".ilist");
    save_carr(ilist,ilist_file);
    ilist_file.close();
    save_content(bwlast,filename + ".bwlast");
    auto occ_file = ofstream(filename + ".occ");
    save_carr(occ,occ_file);
    occ_file.close();
}