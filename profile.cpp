#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <ctime>
#include <algorithm>

struct Gene {
    std::string id;
    std::string name;
    std::string type;
    int start;
    int end;
    std::vector<std::pair<int, int>> exons;
    int count;
    int length;
};

void parseGTF(const std::string &filename, std::vector<std::vector<Gene>> &genelist) {
    std::ifstream file(filename);
    std::string line;
    while (std::getline(file, line)) {
        if (line[0] == '#') continue;
        std::istringstream iss(line); // streaming a line
        std::string chrom, source, feature, start, end, score, strand, frame, attributes;
        iss >> chrom >> source >> feature >> start >> end >> score >> strand >> frame;
        int chrom_p;
        if (chrom == "chrX") {
            chrom_p = 23;
        }
        else if (chrom == "chrY") {
            chrom_p = 24;
        }
        else if (chrom == "chrM") {
            chrom_p = 25;
        }
        else {
            chrom_p = std::stoi(chrom.substr(3));
        }
        int start_p = std::stoi(start);
        int end_p = std::stoi(end);
        if (feature == "gene") {
            std::getline(iss, attributes);
            std::string gene_id, gene_name, gene_type;
            size_t pos;
            if ((pos = attributes.find("gene_id")) != std::string::npos) {
                gene_id = attributes.substr(pos + 9);
                gene_id = gene_id.substr(0, gene_id.find_first_of("\";"));
            }
            if ((pos = attributes.find("gene_name")) != std::string::npos) {
                gene_name = attributes.substr(pos + 11);
                gene_name = gene_name.substr(0, gene_name.find_first_of("\";"));
            }
            if ((pos = attributes.find("gene_type")) != std::string::npos) {
                gene_type = attributes.substr(pos + 11);
                gene_type = gene_type.substr(0, gene_type.find_first_of("\";"));
            }
            Gene newgene = {gene_id, gene_name, gene_type, start_p, end_p, {}, 0, 0};
            genelist[chrom_p].push_back(newgene);
        }
        else if (feature == "exon") {
            int gene_last = genelist[chrom_p].size()-1;
            genelist[chrom_p][gene_last].exons.push_back({start_p, end_p});
            genelist[chrom_p][gene_last].length += (end_p - start_p);
        }
    }
    for (int i=1; i<=25; ++i) {
        std::sort(genelist[i].begin(), genelist[i].end(), [](Gene a, Gene b){return (a.start<b.start) ? 1 : 0;});
        for (int j=0; j<genelist[i].size(); ++j) {
            std::sort(genelist[i][j].exons.begin(), genelist[i][j].exons.end());
        }
    }
    std::cout << '\r' << "[Read Gene] " << "Complete" << "                    \n" << std::flush; 
    // for (Gene g : genelist[1]) {
    //     std::cout << g.id << ' ' << g.name << ' ' << g.type << ' ' << g.start << ' ' << g.end << '\n';
    //     for (auto p : g.exons) {
    //         std::cout << p.first << ' ' << p.second << " / ";
    //     }
    //     std::cout << '\n';
    // }
}

void parseSAM(const std::string &filename, std::vector<std::vector<Gene>> &genelist) {
    std::ifstream file(filename);
    std::string line;
    int read_c = 0;
    while (std::getline(file, line)) {
        if (line[0] == '@') continue;
        std::istringstream iss(line);
        std::string qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual;
        iss >> qname >> flag >> rname >> pos >> mapq >> cigar >> rnext >> pnext >> tlen >> seq >> qual;
        int pos_i = std::stoi(pos);

        // if (pos_i!=6678371)
        //     continue;

        int mapq_i = std::stoi(mapq);
        if (mapq_i < 20)
            continue;
        int chr_i = 0;
        std::string chr_pre = rname.substr(0, rname.find('.'));
        if (chr_pre == "NC_012920")
            chr_i = 25;
        else
            std::istringstream(chr_pre.substr(chr_pre.size()-2)) >> chr_i;
        int gl = -1, gr = genelist[chr_i].size(); // search gene
        while (gr - gl > 1) {
            int gm = (gl + gr) / 2;
            if (pos_i < genelist[chr_i][gm].start)
                gr = gm;
            else
                gl = gm;
        }
        int fg = 0;
        for (int p=gl; p>=0 && p>gl-10; --p) { // search exon
            for (auto ex : genelist[chr_i][p].exons) {
                // std::cout << std::min(pos_i+(int)seq.size(),ex.second) - std::max(pos_i,ex.first) << '\n';
                if ((std::min(pos_i+(int)seq.size(),ex.second)-std::max(pos_i,ex.first))*3 > (int)seq.size()) { // Overlap part > 1/2
                    genelist[chr_i][p].count += 1;
                    // std::cout << chr_i << ' ' << pos_i << ' ' << genelist[chr_i][p].name << '\n';
                    fg = 1;
                    break;
                }
            }
        }
        std::cout <<'\r' << "[Find Gene] " << ++read_c << "                    " << std::flush;
        // if (!fg)
            // std::cout << chr_i << ' ' << pos_i << '\n';
        // std::cout << chr_i << '\n';
    }
    std::cout << '\r' << "[Find Gene] " << "Complete" << "                    \n" << std::flush; 
}

void profile(char** argv) {
    std::time_t start, finish;
    time(&start);
    std::vector<std::vector<Gene>> genelist(30); // genes in each chromosome
    parseGTF(argv[1], genelist);
    parseSAM(argv[2], genelist);
    time(&finish);
    std::cout << '\r' << "[Finish] Total time: " << difftime(finish, start) << " seconds" << "                    \n" << std::flush;
}

int main(int argc, char** argv) {
    profile(argv);
    return 0;
}