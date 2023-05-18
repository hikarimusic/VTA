#include <cstdlib>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>

void grch_nc(char** argv) {
    std::ifstream infile;
    infile.open(argv[1], std::ios::binary);
    infile.seekg(0, std::ios::end);
    std::uint32_t seq_l = infile.tellg();
    char* read(new char[seq_l+1]{});
    infile.seekg(0, std::ios::beg);
    infile.read((char*)read, seq_l);
    infile.close();
    char* write{new char[seq_l+1]{}};
    std::uint32_t read_f{0}, cnt{0};
    for (std::uint32_t i=0; i<=seq_l; ++i) {
        if (!(read[i]))
            break;
        if (read[i]=='>') {
            if (read[i+1]=='N' && read[i+2]=='C')
                read_f = 1;
            else
                read_f = 0;
        }
        if (read_f==1) {
            write[cnt] = read[i];
            cnt += 1;
        }
        if ((i % 1024)==0)
            std::cout << '\r' << std::flush << "[Process] " << i << "/" << seq_l << "                    ";
    }
    std::ofstream outfile;
    outfile.open(std::string(argv[1]).substr(0, 6)+".NC.fa", std::ios::trunc | std::ios::binary);
    outfile.write((char*) write, cnt);
    outfile.close();
    delete[] read;
    delete[] write;
    std::cout << '\r' << std::flush << "[Process] " << "Complete" << "                    " << '\n';
}

int main(int argc, char** argv) {
    grch_nc(argv);
    return 0;
}