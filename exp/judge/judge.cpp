#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

std::string rcseq(std::string s) {
    std::string res;
    for (std::int32_t i=s.size()-1; i>=0; --i) {
        if (s[i]=='A')
            res += 'T';
        else if (s[i]=='C')
            res += 'G';
        else if (s[i]=='G')
            res += 'C';
        else if (s[i]=='T')
            res += 'A';
    }
    return res;
}

void judge() {
    const int len{10000000};
    char* seq{new char[len]{}};
    std::vector<std::string> chr_n{"NC_1", "NC_2", "NC_3", "NC_4", "NC_5", "NC_6", "NC_7", "NC_8", "NC_9", "NC_10"};
    std::ifstream sequence;
    std::ifstream input;
    std::ifstream output;
    std::ifstream submit;
    sequence.open("sequence.fna");
    input.open("input.txt");
    output.open("output.txt");
    submit.open("submit.txt");
    std::string str_i;
    int ni{};
    while (sequence>>str_i) {
        if (str_i[0]!='>') {
            for (int i=0; i<str_i.size(); ++i) {
                seq[ni] = str_i[i];
                ni += 1;
            }
        }
    }
    int qn{};
    int err{};
    std::string qry;
    std::string chr;
    std::int64_t pos{};
    std::string rf;
    std::string aln;
    std::vector<int> aln_i;
    std::vector<char> aln_c;
    int val_m{}, val_p{}, val_d{}, tmp_d{};
    int qs{}, rs{};
    input >> qn;
    for (int qi=0; qi<qn; ++qi) {
        for (int q=0; q<2; ++q) {
            input >> qry;
            submit >> chr >> pos >> rf >> aln;
            if (rf=="R")
                qry = rcseq(qry);
            pos = pos + ((std::find(chr_n.begin(), chr_n.end(), chr)-chr_n.begin())*len/10) - 1;
            aln_i.clear();
            aln_c.clear();
            for (int i=0; i<aln.size(); ++i) {
                if (aln[i]>='A' && aln[i]<='Z') {
                    aln_i.push_back(0);
                    aln_c.push_back(aln[i]);
                }
                else
                    aln_i.back() = aln_i.back() * 10 + (aln[i]-'0');
            }
            val_m = 0;
            val_p = pos;
            val_d = 0;
            qs = 0;
            rs = pos;
            for (int i=0; i<aln_i.size(); ++i) {
                if (aln_c[i]=='M') {
                    for (int j=0; j<aln_i[i]; ++j) {
                        if (seq[rs]!=qry[qs]) {
                            std::cout << "aln: " << i << ' ' << j << ' ' << aln_i[i] << ' ' << aln_c[i] << '\n';
                            std::cout << "M: " << rs-pos << seq[rs] << ' ' << qs << qry[qs] << '\n';
                            val_m = 1;
                            break;
                        }
                        qs += 1;
                        rs += 1;
                    }
                }
                else if (aln_c[i]=='S') {
                    for (int j=0; j<aln_i[i]; ++j) {
                        if (seq[rs]==qry[qs]) {
                            std::cout << "S: " << rs-pos << seq[rs] << ' ' << qs << qry[qs] << '\n';
                            val_m = 1;
                            break;
                        }
                        qs += 1;
                        rs += 1;
                        val_d += 1;
                    }
                }
                else if (aln_c[i]=='D') {
                    for (int j=0; j<aln_i[i]; ++j) {
                        rs += 1;
                        val_d += 1;
                    }
                }
                else {
                    for (int j=0; j<aln_i[i]; ++j) {
                        qs += 1;
                        val_d += 1;
                    }
                }
                if (val_m)
                    break;
            }
            if (val_m) {
                std::cout << "Match error: " << qi << "/" << q << '\n';
                err += 1;
                continue;
            }
            output >> chr >> pos >> rf >> aln;
            pos = pos + ((std::find(chr_n.begin(), chr_n.end(), chr)-chr_n.begin())*len/10) - 1;
            if ((val_p-pos)>10 || (val_p-pos)<-10) {
                std::cout << "Position error: " << qi << "/" << q << '\n';
                err += 1;
                continue;
            }
            aln_i.clear();
            aln_c.clear();
            for (int i=0; i<aln.size(); ++i) {
                if (aln[i]>='A' && aln[i]<='Z') {
                    aln_i.push_back(0);
                    aln_c.push_back(aln[i]);
                }
                else
                    aln_i.back() = aln_i.back() * 10 + (aln[i]-'0');
            }
            tmp_d = 0;
            for (int i=0; i<aln_i.size(); ++i) {
                if (aln_c[i]!='M')
                    tmp_d += aln_i[i];
            }
            if (val_d>tmp_d) {
                std::cout << "Distance error: " << qi << "/" << q << '\n';
                err += 1;
            }
        }
    }
    std::cout << "Error rate: " << err << "/" << qn << '\n';
    sequence.close();
    input.close();
    output.close();
    submit.close();
    delete[] seq;
}

int main() {
    judge();
    return 0;
}