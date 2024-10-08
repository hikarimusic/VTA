#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <ctime>
#include <algorithm>

// Determine the chromosome name, chromosome count, sequence length
void profile(char* seq_f, std::vector<std::string>& chr_n, std::vector<uint32_t>& chr_c, std::uint32_t& len) { // (file, name, count, length)
    std::ifstream infile;
    infile.open(std::string(seq_f)+".chr");
    std::string line;
    while (std::getline(infile, line)) {
        chr_n.push_back(line.substr(0, line.find(' ')));
        std::getline(infile, line);
        chr_c.push_back((std::uint32_t) std::stol(line));
        len += (std::uint32_t) std::stol(line);
    }
    for (int i=1; i<chr_c.size(); ++i) // Prefix sum
        chr_c[i] += chr_c[i-1];
    len = (((len>>6)+1)<<6); // Padding to 64
    infile.close();
}

// Load reference sequence
void loadSeq(char* seq_f, std::uint32_t* seq, std::uint32_t len) { // (file, sequence, length)
    std::cout <<'\r' << "[Load Index] " << "..." << "                    " << std::flush;
    std::ifstream infile;
    infile.open(std::string(seq_f)+".seq", std::ios::binary);
    infile.read((char*) seq, len>>2); // uint32 = 8 bytes
    infile.close();
    std::cout <<'\r' << "[Load Index] " << "Complete" << "                    " << '\n' << std::flush;
}

// Gene object
struct Gene {
    std::string id;
    std::string name;
    std::string type;
    int start;
    int end;
    std::vector<std::pair<int, int>> exons;
    int count;
    int length;
    double tpm;
    double fpkm;
};

// Load gencode annotation
void loadGen(const std::string &filename, std::vector<std::vector<Gene>> &genelist) { // (file, genelist by chromosome)
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
        if (feature == "gene") { // Add new gene
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
        else if (feature == "exon") { // Add new exon
            int gene_last = genelist[chrom_p].size()-1;
            genelist[chrom_p][gene_last].exons.push_back({start_p, end_p});
            genelist[chrom_p][gene_last].length += (end_p - start_p);
        }
    }
    for (int i=1; i<=25; ++i) {
        std::sort(genelist[i].begin(), genelist[i].end(), [](Gene a, Gene b){return (a.start<b.start) ? 1 : 0;}); // Sort genes in each chromosome
        for (int j=0; j<genelist[i].size(); ++j) {
            std::sort(genelist[i][j].exons.begin(), genelist[i][j].exons.end()); // Sort exons in each gene
        }
    }
    std::cout << '\r' << "[Read Gene] " << "Complete" << "                    \n" << std::flush; 
}

// Extract a contiguous segement of nucleotides
std::uint32_t nucs(std::uint32_t* seq, std::uint32_t len, std::uint32_t p, std::uint32_t r) { // (sequence, length, position, number)
    std::uint32_t b = p >> 4; // Find the box of the start nuc
    std::uint32_t h = p & 0b1111; // Find the position of the start nuc (to the current box)
    std::uint32_t t = h + r; // Find the position of the end nuc (to the current box)
    if (p >= len) { // Position is beyond the end of sequence
        return 0;
    }
    else if (p+r >= len) { // Range extends beyond the end of sequence
        t -= 16;
        return (((seq[b]<<(h<<1))>>(h<<1))<<(t<<1)); //  Discard left nucs and make space for remaining dummy As
    }
    else if (t <= 16) { // Range fits within one box
        return ((seq[b]<<(h<<1))>>((16-r)<<1)); // Discard left nucs and align to right
    }
    else { // Range spans two boxes
        t -= 16;
        return  ((((seq[b]<<(h<<1))>>(h<<1))<<(t<<1))+(seq[b+1]>>((16-t)<<1))); // Box1: Discard left nucs and make space for box2. Box2: Discard right nucs
    }
}

void analyze(const std::string &filename, std::uint32_t* seq, std::uint32_t len, std::vector<std::vector<Gene>> &genelist) {
    // std::ifstream file(filename);
    // std::string line;
    // int read_c = 0;
    // std::string cur_q = ""; // current query, for checking paired-reads
    // while (std::getline(file, line)) {
    //     if (line[0] == '@') continue;
    //     std::istringstream iss(line);
    //     std::string qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual;
    //     iss >> qname >> flag >> rname >> pos >> mapq >> cigar >> rnext >> pnext >> tlen >> seq >> qual;
    //     int pos_i = std::stoi(pos);

    //     // if (pos_i!=6678371)
    //     //     continue;

    //     int mapq_i = std::stoi(mapq);
    //     if (mapq_i < 20)
    //         continue;
    //     int chr_i = 0;
    //     std::string chr_pre = rname.substr(0, rname.find('.'));
    //     if (chr_pre == "NC_012920")
    //         chr_i = 25;
    //     else
    //         std::istringstream(chr_pre.substr(chr_pre.size()-2)) >> chr_i;
    //     int gl = -1, gr = genelist[chr_i].size(); // search gene
    //     while (gr - gl > 1) {
    //         int gm = (gl + gr) / 2;
    //         if (pos_i < genelist[chr_i][gm].start)
    //             gr = gm;
    //         else
    //             gl = gm;
    //     }
    //     // int fg = 0;
    //     int find_fg = 0; // match 1 gene, skip other genes
    //     // std::cout << chr_i << ' ' << pos_i << '\n';
    //     for (int p=gl; p>=0 && p>gl-10; --p) { // search exon
    //         if (pos_i > genelist[chr_i][p].end)
    //             continue;
    //         for (auto ex : genelist[chr_i][p].exons) {
    //             // std::cout << std::min(pos_i+(int)seq.size(),ex.second) - std::max(pos_i,ex.first) << '\n';
    //             if ((std::min(pos_i+static_cast<int>(seq.size()),ex.second)-std::max(pos_i,ex.first))*3 > static_cast<int>(seq.size())) { // Overlap part > 1/2
    //                 if (qname == cur_q && qname != "") {
    //                     genelist[chr_i][p].count += 1;
    //                     cur_q = "";
    //                     // std::cout << genelist[chr_i][p].name << '\n';
    //                 }
    //                 else {
    //                     cur_q = qname;
    //                 }
    //                 // genelist[chr_i][p].count += 1;
    //                 // std::cout << chr_i << ' ' << pos_i << ' ' << genelist[chr_i][p].name << '\n';
    //                 // fg = 1;
    //                 find_fg = 1;
    //                 break;
    //             }
    //         }
    //         if (find_fg)
    //             break;
    //     }
    //     std::cout <<'\r' << "[Find Gene] " << ++read_c << "                    " << std::flush;

    //     // if (!fg)
    //         // std::cout << chr_i << ' ' << pos_i << '\n';
    //     // std::cout << chr_i << '\n';
    // }
    std::cout << '\r' << "[Find Gene] " << "Complete" << "                    \n" << std::flush; 
}

// void save(const std::string &filename, std::vector<std::vector<Gene>> &genelist) {
//     std::ofstream file(filename);
//     int frag_sum = 0;
//     double rpk_sum = 0;
//     for (int i=1; i<=25; ++i) { // total fragments, total reads per kilobase
//         for (Gene &g : genelist[i]) {
//             frag_sum += g.count;
//             rpk_sum += static_cast<double>(g.count) / (static_cast<double>(g.length) / 1000.0);
//         }
//     }
//     for (int i=1; i<=25; ++i) { // tpm, rpkm
//         for (Gene &g : genelist[i]) {
//             double rpk = static_cast<double>(g.count) / (static_cast<double>(g.length) / 1000.0);
//             g.tpm = rpk / (rpk_sum / 1000000.0);
//             g.fpkm = rpk / (frag_sum / 1000000.0);
//         }
//     }
//     file << "gene_id\tgene_name\tgene_type\tgene_count\ttpm\tfpkm\n";
//     for (int i=1; i<=25; ++i) {
//         for (Gene &g : genelist[i]) {
//             file << g.id << "\t" << g.name << "\t" << g.type << "\t" << g.count << "\t" << g.tpm << "\t" << g.fpkm << "\n";
//         }
//     }

//     std::vector<std::pair<double,std::string>> gex;
//     for (int i=1; i<=25; ++i) { // test
//         for (Gene &g : genelist[i]) {
//             gex.push_back({g.tpm, g.name});
//         }
//     }
//     sort(gex.begin(), gex.end(), std::greater<std::pair<double, std::string>>());
//     for (int i=0; i<100; ++i) {
//         std::cout << gex[i].se  cond << '\t' << gex[i].first << '\n';
//     }

//     file.close();
//     std::cout << '\r' << "[Save Profile] " << "Complete" << "                    \n" << std::flush;
// }

void variant(char** argv) {
    std::time_t start, finish;
    time(&start);
    std::vector<std::string> chr_n;
    std::vector<std::uint32_t> chr_c;
    std::uint32_t len{};
    profile(argv[1], chr_n, chr_c, len);
    std::uint32_t* seq{new std::uint32_t[len>>4]{}};
    loadSeq(argv[1], seq, len);
    std::vector<std::vector<Gene>> genelist(30);
    loadGen(argv[2], genelist);

    analyze(argv[3], seq, len, genelist);
    // save(argv[4], genelist);
    
    time(&finish);
    std::cout << '\r' << "[Finish] Total time: " << difftime(finish, start) << " seconds" << "                    \n" << std::flush;
}

int main(int argc, char** argv) {
    variant(argv);
    return 0;
}