#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

void problem() {
    std::cout << '\r' << "[Sequence] " << "Generating" << "                    " << std::flush;
    const int len{10000000};
    const int qn{1000000};
    const int ql{100};
    char* seq{new char[len]{}};
    int rand{};
    char itc[8]{'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T'};
    std::map<char, int> cti = {{'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}};
    for (int i=0; i<len; ++i)
        seq[i] = itc[std::rand()%4];
    std::ofstream seq_f;
    seq_f.open("sequence.fna", std::ios::trunc);
    for (int i=0; i<10; ++i) {
        seq_f << ">NC_" << i+1 << '\n';
        for (int j=0; j<len/1000; ++j) {
            seq_f.write((char*)seq+(len/10)*i+100*j, 100);
            seq_f << '\n';
        }
    }
    seq_f.close();
    std::cout << '\r' << "[Sequence] " << "Complete" << "                    \n" << std::flush;
    std::string qry[2];
    int pos[2]{};
    std::vector<int> aln_i[2];
    std::vector<char> aln_c[2];
    std::ofstream input;
    std::ofstream output;
    input.open("input.txt", std::ios::trunc);
    output.open("output.txt", std::ios::trunc);
    input << qn << '\n';
    int insl{}, head{}, tail{}, qs{}, rs{}, rf[2]{};
    for (int qi=0; qi<qn; ++qi) {
        insl = (std::rand() % 600) + 201;
        head = std::rand() % (len-insl+1);
        tail = head + insl - 1;
        pos[0] = head;
        pos[1] = tail;
        for (int q=0; q<2; ++q) {
            qry[q].clear();
            aln_i[q].clear();
            aln_c[q].clear();
            qs = 0;
            rs = pos[q];
            while (qs<ql) {
                rand = std::rand() % 1000;
                if (rand<10) {
                    qry[q].push_back(itc[(q?(3-cti[seq[rs]]):cti[seq[rs]])+std::rand()%3+1]);
                    if (aln_c[q].empty() || aln_c[q].back()!='S') {
                        aln_i[q].push_back(1);
                        aln_c[q].push_back('S');
                    }
                    else
                        aln_i[q].back() += 1;
                    qs += 1;
                    rs += (q)?(-1):1;
                }
                else if (rand<15) {
                    if (aln_c[q].empty() || aln_c[q].back()!='D') {
                        aln_i[q].push_back(1);
                        aln_c[q].push_back('D');
                    }
                    else    
                        aln_i[q].back() += 1;
                    rs += (q)?(-1):1;
                }
                else if (rand<20) {
                    qry[q].push_back(itc[std::rand()%4]);
                    if (aln_c[q].empty() || aln_c[q].back()!='I') {
                        aln_i[q].push_back(1);
                        aln_c[q].push_back('I');
                    }
                    else
                        aln_i[q].back() += 1;
                    qs += 1;
                }
                else {
                    qry[q].push_back(q?itc[3-cti[seq[rs]]]:seq[rs]);
                    if (aln_c[q].empty() || aln_c[q].back()!='M') {
                        aln_i[q].push_back(1);
                        aln_c[q].push_back('M');
                    }
                    else    
                        aln_i[q].back() += 1;
                    qs += 1;
                    rs += (q)?(-1):1;
                }
            }
            if (q) {
                pos[q] = rs + 1;
                std::reverse(aln_i[q].begin(), aln_i[q].end());
                std::reverse(aln_c[q].begin(), aln_c[q].end());
            }
        }
        if (std::rand()&1) {
            rf[0] = 0;
            rf[1] = 1;
        }
        else {
            rf[0] = 1;
            rf[1] = 0;
        }
        for (int q=0; q<2; ++q) {
            input << qry[rf[q]] << ' ';
            input << (q?'\n':' ');
            output << "NC_" << pos[rf[q]]/(len/10)+1 << ' ' << pos[rf[q]]%(len/10)+1 << ' ';
            output << (rf[q]?'R':'F') << ' ';
            for (int i=0; i<aln_i[rf[q]].size(); ++i)
                output << aln_c[rf[q]][i] << aln_i[rf[q]][i];
            output << (q?'\n':' ');
        }
        if (qi%1024==0)
            std::cout << '\r' << "[Query] " << qi << '/' << qn << "                    " << std::flush;
    }
    input.close();
    output.close();
    std::cout << '\r' << "[Query] " << "Complete" << "                    \n" << std::flush;
    delete[] seq;
}

int main() {
    problem();
    return 0;
}