#include <bits/stdc++.h> //reverse
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <math.h>       /* log2 */


class Read {
    public:
        int index;
        int pos;
        std::string seq;

        Read(): pos(-1), seq("") { }
        Read(int _index, int _pos, std::string _seq): index(_index), pos(_pos), seq(_seq) { }
        Read & operator=(const Read & other) {
            pos = other.pos;
            seq = other.seq;
            return *this;
        }
};

Read empty_read = Read();
long long int nb_comparison = 0;
long long int estimated_nb_comparison = 0;
int last_percent = -1;

class Suffix {
    public:
        bool is_in_ref;
        int pos_in_seq;
        std::string & ref;
        Read * read;

        Suffix(bool _in_ref,
                int _pos_in_seq,
                std::string & _ref):
            is_in_ref(_in_ref),
            pos_in_seq(_pos_in_seq),
            ref(_ref),
            read(&empty_read){ }

        Suffix(bool _in_ref,
                int _pos_in_seq,
                std::string & _ref,
                Read & _read):
            is_in_ref(_in_ref),
            pos_in_seq(_pos_in_seq),
            ref(_ref),
            read(&_read) { }

        Suffix & operator=(const Suffix & other) {
        is_in_ref = other.is_in_ref;
        pos_in_seq = other.pos_in_seq;
        ref = other.ref;
        read = other.read;
        return *this;
        }

        char get_char(int i) const {
            if (is_in_ref) {
                if (i+pos_in_seq < ref.size())
                    return ref[i+pos_in_seq];
                else
                    return '$';
            }
            else {
                if (i+pos_in_seq < read->seq.size())
                    return read->seq[i+pos_in_seq];
                else {
                    int i_ref = read->pos + i + pos_in_seq - read->seq.size();
                    if (i_ref < ref.size())
                        return ref[i_ref];
                    else
                        return '$';
                }
            }
        }

        bool operator<(const Suffix & other) const {
            nb_comparison++;
            if (int(((float)nb_comparison / estimated_nb_comparison)*100) != last_percent ){
                std::cout << int(((float) nb_comparison / estimated_nb_comparison)*100) << "%" << std::endl;
                last_percent = int(((float)nb_comparison / estimated_nb_comparison)*100);
            }
            int i = 0;
            while (true) {
                char c = get_char(i);
                char other_c = other.get_char(i);
                if (c == '$' && other_c == '$'){
                    if (is_in_ref && other.is_in_ref)
                        return pos_in_seq < other.pos_in_seq;
                    else if (!is_in_ref && !other.is_in_ref)
                        return read->index < other.read->index;
                    else if (is_in_ref && !other.is_in_ref)
                        return true;
                    else
                        return false;
                }
                if (c == '$')
                    return true;
                if (other_c == '$')
                    return false;
                if (c != other_c)
                    return c < other_c;
                i++;
            }
        }

};

std::ostream &operator<<(std::ostream &os, Read const &r) {
    return os << r.pos << " " << r.seq;
}

std::ostream &operator<<(std::ostream &os, Suffix const &s) {
    if (s.is_in_ref){
        return os << "ref " << s.pos_in_seq;
    }
    else {
        return os << "read " << s.pos_in_seq << " " << *s.read;
    }
}

int main(int argc,char* argv[]) {
    if (argc==1) {
        std::cout << "You must give the input file as an argument after the program name." << std::endl;
        return 1;
    }
    else {
        int w = 10;
        std::string ref;
        std::vector<Read> reads;
        std::ifstream inFile;
        inFile.open(argv[1]);
        if (!inFile) {
            std::cout << "Unable to open file";
            exit(1); // terminate with error
        }

        //parse custom format
        inFile >> w;
        inFile >> ref;
        std::reverse(ref.begin(),ref.end());
        int pos = 0;
        std::string seq = "";
        int total_size_read = 0;
        std::string line;
        int index = 0;
        while (std::getline(inFile,line)){
            std::stringstream line_stream(line);
            if (line_stream >> pos >> seq) {
            reverse(seq.begin(),seq.end());
            if (pos == 0) reads.push_back(Read(index,ref.size()-pos,seq));
            else reads.push_back(Read(index,ref.size()-pos-w,seq));
            index++;
            total_size_read = total_size_read + seq.size();
            }
        }
        std::cout << "Finished parsing input" << std::endl;
        std::cout << "sum: " <<  ref.size() + total_size_read << std::endl;

        //creation of the sa
        std::vector<Suffix> sa;
        std::vector<std::vector<int>> start_of_read(ref.size()+1);
        for(int i = 0; i < ref.size(); i++) {
            sa.push_back(Suffix(true,i,ref));
        }
        for(int r = 0; r < reads.size(); r++) {
            start_of_read[reads[r].pos].push_back(r);
            for(int i= 0; i < reads[r].seq.size(); i++){
                sa.push_back(Suffix(false,i,ref,reads[r]));
            }
        }
        /*
        std::cout << "Finished adding to sa" << std::endl;
        std::cout << "size ref: " << ref.size() << std::endl;
        std::cout << "total size of reads: " << total_size_read << std::endl;
        std::cout << "size sa: " << sa.size() << std::endl;
        */

        estimated_nb_comparison = sa.size()*log2(sa.size());
        std::sort(sa.begin(),sa.end());
        std::cout << "Finished ordering sa" << std::endl;
        std::cout << "sa size: " << sa.size() << std::endl;
        std::cout << "number of reads: " << reads.size() << std::endl;

        //creation of the bwt
        std::vector<char> bwt;
        bwt.push_back(ref.back());
        for(int i = 0; i < start_of_read[ref.size()].size(); i++)
            bwt.push_back(reads[i].seq.back());
        for(int i = 0; i < sa.size(); i++) {
            /*if (i== 838 || i == 839){
                std::cout << "is_in_ref: " << sa[i].is_in_ref << std::endl;
                std::cout << "pos in read: "<< sa[i].pos_in_seq << std::endl;
                std::cout << "index: " << sa[i].read->index << std::endl;
                std::cout << "pos in reference: " << sa[i].read->pos << std::endl;
                std::cout << "index: " << sa[i].read->seq << std::endl;
            }*/
            if (sa[i].pos_in_seq > 0){
                bwt.push_back(sa[i].get_char(-1));
                if (sa[i].is_in_ref){
                for(int j = 0; j < start_of_read[sa[i].pos_in_seq].size(); j++){
                    bwt.push_back(reads[start_of_read[sa[i].pos_in_seq][j]].seq.back());
                }
                }
            }
        }
        std::cout << "Finished building bwt" << std::endl;
        std::cout << "bwt size: " << bwt.size() << std::endl;

        //output the result
        //std::cout <<"BWT:" << std::endl;
        std::ofstream bwt_file;
        bwt_file.open((std::string) argv[1] + ".bwt");
        for(int i = 0; i < bwt.size(); i++)
            bwt_file << bwt[i];
        bwt_file.close();
        char last = bwt[0];
        int nb_last = 1;
        int nb_runs = 1;
        for(int i = 1; i < bwt.size(); i++) {
            if (bwt[i] == last)
                nb_last++;
            else {
                //std::cout << last << nb_last;
                last = bwt[i];
                nb_last = 1;
                nb_runs++;
            }
        }
        //std::cout << std::endl;
        std::cout << "Number of runs: " << nb_runs << std::endl;
        std::cout << "average size of run: " << sa.size()/nb_runs << std::endl;

    }
    return 0;
}
