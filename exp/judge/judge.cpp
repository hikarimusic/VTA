#include <algorithm>
#include <cstdint>
#include <cstdlib>
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

bool parse(std::ifstream& infile, std::vector<std::string>& line) {
    std::string str_i;
    if (std::getline(infile, str_i)) {
        line.clear();
        line.push_back("");
        for (int i=0; i<str_i.size(); ++i) {
            if (str_i[i]=='\t')
                line.push_back("");
            else
                line.back().push_back(str_i[i]);
        }
        return true;
    }
    else
        return false;
}

void judge(char** argv) {
    std::ifstream sequence;
    std::ifstream output;
    std::ifstream submit;
    sequence.open("sequence.fna");
    output.open("output.sam");
    submit.open(std::string(argv[1]));
    sequence.seekg(0, std::ios::end);
    std::uint32_t len = sequence.tellg();
    char* seq{new char[len]{}};
    sequence.seekg(0, std::ios::beg);
    std::vector<std::string> chr_n;
    std::vector<std::uint32_t> chr_c;
    std::string str_i;
    std::uint32_t ni{};
    while (sequence>>str_i) {
        if (str_i[0]=='>') {
            chr_n.push_back(str_i.substr(1, str_i.size()-1));
            if (ni>0)
                chr_c.push_back(ni);
        }
        else {
            for (int i=0; i<str_i.size(); ++i) {
                seq[ni] = str_i[i];
                ni += 1;
            }
        }
    }
    chr_c.push_back(ni);
    std::cout << '\r' << "[Load] " << "Complete" << "                    \n" << std::flush;
    int qn{};
    int err[2]{};
    std::vector<std::string> line[2];
    std::string qry_tp{""};
    int cig_tp[2]{};
    while (parse(output, line[0])) {
        if (line[0][0][0]!='@')
            break;
    }
    while (parse(submit, line[1])) {
        if (line[1][0][0]!='@')
            break;
    }
    // std::cout << line[0][9] << '\n' << line[1][9] << '\n';
    do {
        qn += 1;
        // std::cout << line[0][9] << '\n' << line[1][9] << '\n';
        // std::cout << cig_tp[0] << ' ' << cig_tp[1] << '\n';
        // std::cout << line[1][1] << ' ' << ((std::stoi(lin[1][1])&128)/128) << '\n';
        while ((qry_tp==line[1][0]) && (cig_tp[(std::stoi(line[1][1])&128)/128])) {
            std::cout << '\r' << "    E. Duplicated: " << line[1][0] << "                    \n" << std::flush;
            parse(submit, line[1]);
        }
        if (qry_tp!=line[1][0]) {
            qry_tp = line[1][0];
            cig_tp[0] = 0;
            cig_tp[1] = 0;
        }
        cig_tp[(std::stoi(line[1][1])&128)/128] = 1;
        while ((qry_tp!=line[0][0]) || (std::stoi(line[0][1])&128!=std::stoi(line[1][1])&128)) {
            std::cout << '\r' << "    E. Unmaped: " << line[0][0] << "                    \n" << std::flush;
            err[0] += 1;
            parse(output, line[0]);
            qn += 1;
        }
        if (line[0][2]!=line[1][2] || std::abs(std::stol(line[0][3])-std::stol(line[1][3]))>100) {
            std::cout << '\r' << "    E. Position: " << line[0][0] << "                    \n" << std::flush;
            err[0] += 1;
            continue;
        }
        std::string qry;
        std::string chr;
        std::int64_t pos{};
        std::string aln;
        std::vector<int> aln_i;
        std::vector<char> aln_c;
        int edit[2]{};
        for (int q=0; q<2; ++q) {
            edit[q] = 0;
            qry = line[q][9];
            chr = line[q][2];
            pos = stol(line[q][3]);
            aln = line[q][5];
            aln_i.clear();
            aln_c.clear();
            aln_i.push_back(0);
            for (int i=0; i<aln.size(); ++i) {
                if (aln[i]>='A' && aln[i]<='Z') {
                    aln_c.push_back(aln[i]);
                    aln_i.push_back(0);
                }
                else
                    aln_i.back() = aln_i.back() * 10 + (aln[i]-'0');
            }
            aln_i.pop_back();
            int r = std::find(chr_n.begin(), chr_n.end(), chr)-chr_n.begin();
            if (r<chr_n.size())
                pos = pos + (r?chr_c[r-1]:0) - 1;
            else
                pos -= 1;
            int qs = 0;
            int rs = pos;
            for (int i=0; i<aln_c.size(); ++i) {
                if (aln_c[i]=='M') {
                    for (int j=0; j<aln_i[i]; ++j) {
                        if (seq[rs]!=qry[qs])
                            edit[q] += 1;
                        qs += 1;
                        rs += 1;
                    }
                }
                else if (aln_c[i]=='I' || aln_c[i]=='S') {
                    for (int j=0; j<aln_i[i]; ++j) {
                        edit[q] += 1;
                        qs += 1;
                    }
                }
                else if (aln_c[i]=='D' || aln_c[i]=='N') {
                    for (int j=0; j<aln_i[i]; ++j) {
                        edit[q] += 1;
                        rs += 1;
                    }
                }
                else {
                    for (int j=0; j<aln_i[i]; ++j) {
                        edit[q] += 1;
                        qs += 1;
                        rs += 1;
                    }
                }
            }

        }
        if (edit[1]>edit[0]) {
            std::cout << '\r' << "    E. Distance: " << line[0][0] << "                    \n" << std::flush;
            err[1] += 1;
            continue;
        }
        std::cout << '\r' << "[Judge] " << qn << "                    " << std::flush;
    } while (parse(output, line[0]) && parse(submit, line[1]));

    // while (std::getline(submit, str_i)) {
    //     parse(str_i, line[0]);
    //     std::getline(output, str_i);
    //     parse(str_i, line[1]);
    //     while (line)
    //     std::string qry = line[9];
    //     std::string chr = line[2];
    //     std::int64_t pos = std::stol(line[3]);
    //     int r = std::find(chr_n.begin(), chr_n.end(), chr)-chr_n.begin();
    //     if (r<chr_n.size())
    //         pos = pos + (r?chr_c[r-1]:0) - 1;
    // }

    // int qn{};
    // int err{};
    // std::string qry;
    // std::string chr;
    // std::int64_t pos{};
    // std::string rf;
    // std::string aln;
    // std::vector<int> aln_i;
    // std::vector<char> aln_c;
    // int val_m{}, val_p{}, val_d{}, tmp_d{};
    // int qs{}, rs{};
    // input >> qn;
    // for (int qi=0; qi<qn; ++qi) {
    //     for (int q=0; q<2; ++q) {
    //         input >> qry;
    //         submit >> chr >> pos >> rf >> aln;
    //         if (rf=="R")
    //             qry = rcseq(qry);
    //         pos = pos + ((std::find(chr_n.begin(), chr_n.end(), chr)-chr_n.begin())*len/10) - 1;
    //         aln_i.clear();
    //         aln_c.clear();
    //         for (int i=0; i<aln.size(); ++i) {
    //             if (aln[i]>='A' && aln[i]<='Z') {
    //                 aln_i.push_back(0);
    //                 aln_c.push_back(aln[i]);
    //             }
    //             else
    //                 aln_i.back() = aln_i.back() * 10 + (aln[i]-'0');
    //         }
    //         val_m = 0;
    //         val_p = pos;
    //         val_d = 0;
    //         qs = 0;
    //         rs = pos;
    //         for (int i=0; i<aln_i.size(); ++i) {
    //             if (aln_c[i]=='M') {
    //                 for (int j=0; j<aln_i[i]; ++j) {
    //                     if (seq[rs]!=qry[qs]) {
    //                         std::cout << "aln: " << i << ' ' << j << ' ' << aln_i[i] << ' ' << aln_c[i] << '\n';
    //                         std::cout << "M: " << rs-pos << seq[rs] << ' ' << qs << qry[qs] << '\n';
    //                         val_m = 1;
    //                         break;
    //                     }
    //                     qs += 1;
    //                     rs += 1;
    //                 }
    //             }
    //             else if (aln_c[i]=='S') {
    //                 for (int j=0; j<aln_i[i]; ++j) {
    //                     if (seq[rs]==qry[qs]) {
    //                         std::cout << "S: " << rs-pos << seq[rs] << ' ' << qs << qry[qs] << '\n';
    //                         val_m = 1;
    //                         break;
    //                     }
    //                     qs += 1;
    //                     rs += 1;
    //                     val_d += 1;
    //                 }
    //             }
    //             else if (aln_c[i]=='D') {
    //                 for (int j=0; j<aln_i[i]; ++j) {
    //                     rs += 1;
    //                     val_d += 1;
    //                 }
    //             }
    //             else {
    //                 for (int j=0; j<aln_i[i]; ++j) {
    //                     qs += 1;
    //                     val_d += 1;
    //                 }
    //             }
    //             if (val_m)
    //                 break;
    //         }
    //         if (val_m) {
    //             std::cout << "Match error: " << qi << "/" << q << '\n';
    //             err += 1;
    //             continue;
    //         }
    //         output >> chr >> pos >> rf >> aln;
    //         pos = pos + ((std::find(chr_n.begin(), chr_n.end(), chr)-chr_n.begin())*len/10) - 1;
    //         if ((val_p-pos)>10 || (val_p-pos)<-10) {
    //             std::cout << "Position error: " << qi << "/" << q << '\n';
    //             err += 1;
    //             continue;
    //         }
    //         aln_i.clear();
    //         aln_c.clear();
    //         for (int i=0; i<aln.size(); ++i) {
    //             if (aln[i]>='A' && aln[i]<='Z') {
    //                 aln_i.push_back(0);
    //                 aln_c.push_back(aln[i]);
    //             }
    //             else
    //                 aln_i.back() = aln_i.back() * 10 + (aln[i]-'0');
    //         }
    //         tmp_d = 0;
    //         for (int i=0; i<aln_i.size(); ++i) {
    //             if (aln_c[i]!='M')
    //                 tmp_d += aln_i[i];
    //         }
    //         if (val_d>tmp_d) {
    //             std::cout << "Distance error: " << qi << "/" << q << '\n';
    //             err += 1;
    //         }
    //     }
    // }
    // std::cout << "Error rate: " << err << "/" << qn << '\n';
    std::cout << '\r' << "[Judge] " << "Complete" << "                    \n" << std::flush; 
    std::cout << '\r' << "[Result] " << "                    \n" << std::flush;
    std::cout << '\r' << "    Mapping error: " << err[0] << '/' << qn << '\n';
    std::cout << '\r' << "    Alignment error: " << err[1] << '/' << qn << '\n';
    sequence.close();
    output.close();
    submit.close();
    delete[] seq;
}

int main(int argc, char** argv) {
    judge(argv);
    return 0;
}