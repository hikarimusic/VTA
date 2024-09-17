#include <cstdlib>
#include <cstdint>
#include <ctime>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

void grch_nc(char* seq_f) {
    std::ifstream infile;
    infile.open(seq_f, std::ios::binary);
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
            std::cout << '\r' << std::flush << "[Genome] " << i << "/" << seq_l << "                    ";
    }
    std::ofstream outfile;
    outfile.open(seq_f, std::ios::trunc | std::ios::binary);
    outfile.write((char*) write, cnt);
    outfile.close();
    delete[] read;
    delete[] write;
    std::cout << '\r' << std::flush << "[Genome] " << "Complete" << "                    " << '\n';
}

void gencode_g(char* gen_f) {
    std::cout << gen_f << '\n';
    std::ifstream infile;
    infile.open(std::string(gen_f));

    if (!infile.is_open()) {
        std::cerr << "Error: Cannot open the file " << gen_f << '\n';
        return;
    }

    std::string line;
    std::vector<std::string> gene;
    int i = 0;
    std::cout << "hello\n";
    while (std::getline(infile, line)) {
        // Use a string stream to split the line
        std::istringstream stream(line);
        std::string field;
        int column_index = 0;
        bool is_gene = false;

        // Split the line by tabs
        while (std::getline(stream, field, '\t')) {
            // Check the third element (index 2)
            if (column_index == 2 && field == "gene") {
                is_gene = true;
                break; // No need to continue parsing
            }
            column_index++;
        }

        // If the third element is "gene", push the line into the vector
        if (is_gene) {
            gene.push_back(line);
        }

        i += 1;
        if ((i % 1024) == 0) {
            std::cout << '\r' << std::flush << "[Gencode] " << i;
        }
    }
    infile.close();
    std::ofstream outfile;
    outfile.open(gen_f, std::ios::trunc);
    for (std::string s : gene) {
        outfile << s << '\n';
    }
    outfile.close();
    std::cout << '\r' << std::flush << "[Gencode] " << "Complete" << "                    " << '\n';
    std::cout << gene.size() << '\n';
}

void clean(char** argv) {
    std::time_t start, finish;
    time(&start);
    // grch_nc(argv[1]);
    gencode_g(argv[2]);
    time(&finish);
    std::cout << '\r' << "[Finish] Total time: " << difftime(finish, start) << " seconds" << "                    \n" << std::flush;
}

int main(int argc, char** argv) {
    clean(argv);
    return 0;
}