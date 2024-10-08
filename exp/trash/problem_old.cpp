#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
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

void locate(std::string& chr, std::int64_t& pos, std::vector<std::string> chr_n, std::vector<std::int64_t> chr_c) {
    int l{}, m{}, r{};
    r = (int) (chr_n.size()-1);
    while (r>l) {
        m = (l + r) / 2;
        if (pos<chr_c[m])
            r = m;
        else
            l = m + 1;
    }
    chr = chr_n[r];
    pos = pos - (r?chr_c[r-1]:0) + 1;
}

void problem(char** argv) {
    std::int64_t len{4294967296};
    const std::int64_t qn{1000000};
    const std::int64_t ql{100};
    const std::int64_t var_s{10};
    const std::int64_t var_d{5};
    const std::int64_t var_i{5};
    char* seq{new char[len]{}};
    std::vector<std::string> chr_n;
    std::vector<std::int64_t> chr_c;
    std::ifstream seq_f;
    seq_f.open(argv[1]);
    std::string str_i;
    std::int64_t ni{};
    while (seq_f>>str_i) {
        if (str_i[0]=='>') {
            chr_n.push_back(str_i.substr(1, str_i.size()-1));
            if (ni>0)
                chr_c.push_back(ni);
        }
        else {
            for (std::int64_t i=0; i<str_i.size(); ++i) {
                seq[ni] = (str_i[i]>'Z')?(str_i[i]-32):str_i[i];
                ni += 1;
            }
        }
        if ((ni % 1024)==0)
            std::cout << '\r' << "[Sequence] " << ni << "                    " << std::flush;
    }
    chr_c.push_back(ni);
    len = ni;
    seq_f.close();
    std::cout << '\r' << "[Sequence] " << "Complete" << "                    \n" << std::flush;
    std::ofstream input[2];
    std::ofstream output;
    input[0].open("input1.fq", std::ios::trunc);
    input[1].open("input2.fq", std::ios::trunc);
    output.open("output.sam", std::ios::trunc);
    for (int i=0; i<chr_n.size(); ++i) {
        output << "@SQ" << '\t';
        output << "SN:" << chr_n[i] << '\t';
        output << "LN:" << (i?(chr_c[i]-chr_c[i-1]):chr_c[i]) << '\n';
    }
    output << "@PG" << '\t' << "ID:NA" << '\t' << "PN:NA" << '\t' << "VN:NA" << '\t';
    output << "CL:" << "./problem" << '\n';
    std::string qry[2];
    std::string chr[2];
    std::int64_t pos[2]{};
    std::vector<int> aln_i[2];
    std::vector<char> aln_c[2];
    std::int64_t insl{}, head{}, tail{}, qs{}, rs{}, rf[2]{}, rand{};
    char itc[8]{'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T'};
    std::map<char, int> cti = {{'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}};
    for (int qi=0; qi<qn; ++qi) {
        do {
            insl = (std::rand() % 600) + 201;
            head = (((std::int64_t)std::rand())*((std::int64_t)RAND_MAX+1)+((std::int64_t)std::rand())) % (len-insl+1);
            tail = head + insl - 1;
            pos[0] = insl;
            pos[1] = tail;
            locate(chr[0], pos[0], chr_n, chr_c);
            locate(chr[1], pos[1], chr_n, chr_c);
        } while (chr[0]!=chr[1]);
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
                if (rand<var_s) {
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
                else if (rand<var_s+var_d) {
                    if (aln_c[q].empty() || aln_c[q].back()!='D') {
                        aln_i[q].push_back(1);
                        aln_c[q].push_back('D');
                    }
                    else    
                        aln_i[q].back() += 1;
                    rs += (q)?(-1):1;
                }
                else if (rand<var_s+var_d+var_i) {
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
            input[q] << "@SIM" << std::string(7-std::to_string(qi+1).size(), '0')+std::to_string(qi+1) << ' ' << q+1 << '\n';
            input[q] << qry[rf[q]] << '\n' << '+' << '\n' << std::string(ql, 'E') << '\n';
            output << "SIM" << std::string(7-std::to_string(qi+1).size(), '0')+std::to_string(qi+1) << '\t';
            if (rf[0])
                output << (q?163:83) << '\t';
            else
                output << (q?147:99) << '\t';
            locate(chr[rf[q]], pos[rf[q]], chr_n, chr_c);
            output << chr[rf[q]] << '\t';
            output << pos[rf[q]] << '\t';
            output << 255 << '\t';
            std::vector<int> nla_i;
            std::vector<char> nla_c;
            for (int i=0; i<aln_i[rf[q]].size(); ++i) {
                if (aln_c[rf[q]][i]=='M' || aln_c[rf[q]][i]=='S') {
                    if (i==0 || (aln_c[rf[q]][i-1]!='M' && aln_c[rf[q]][i-1]!='S')) {
                        nla_i.push_back(aln_i[rf[q]][i]);
                        nla_c.push_back('M');
                    }
                    else {
                        nla_i.back() += aln_i[rf[q]][i];
                    }
                }
                else {
                    nla_i.push_back(aln_i[rf[q]][i]);
                    nla_c.push_back(aln_c[rf[q]][i]);
                }
            }
            for (int i=0; i<nla_i.size(); ++i)
                output << nla_i[i] << nla_c[i];
            output << '\t';
            output << "=" << '\t';
            output << pos[rf[1-q]]%(len/10)+1 << '\t';
            std::int64_t tmpl{};
            tmpl = pos[1] - pos[0];
            for (int i=0; i<aln_i[1].size(); ++i) {
                if (aln_c[1][i]!='I')
                    tmpl += aln_i[1][i];
            }
            if (rf[q])
                tmpl *= -1;
            output << tmpl << '\t';
            output << (rf[q]?rcseq(qry[rf[q]]):qry[rf[q]]) << '\t';
            output << std::string(ql, 'E') << '\t';
            output << "aln:";
            for (int i=0; i<aln_i[rf[q]].size(); ++i)
                output << aln_c[rf[q]][i] << aln_i[rf[q]][i];
            output << '\n';
        }
        if (qi%1024==0)
            std::cout << '\r' << "[Query] " << qi << '/' << qn << "                    " << std::flush;
    }
    input[0].close();
    input[1].close();
    output.close();
    std::cout << '\r' << "[Query] " << "Complete" << "                    \n" << std::flush;
    delete[] seq;
}

int main(int argc, char** argv) {
    problem(argv);
    return 0;
}