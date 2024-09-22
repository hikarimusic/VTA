#include <cstdlib>
#include <cstdint>
#include <ctime>
#include <fstream>
#include <iostream>
#include <string>

// Determine length of reference sequence
void profile(char* seq_f, std::uint32_t& len) { // (file, length)
    std::ifstream infile;
    infile.open(seq_f, std::ios::binary);
    infile.seekg(0, std::ios::end);
    len = infile.tellg();
    len += 1024;
    infile.close();
}

// Encode nucleotide into packed form
void write_acgt(std::uint32_t* seq, std::uint32_t pos, char nuc) { // (sequence, position, nucleotide)
    std::uint32_t box = pos >> 4; // Find the box of nucleotide (1 nuc = 2 bit, so 16 nucs = uint32)
    std::uint32_t pnt = (0b1111 - pos & 0b1111) << 1; // Find the position in box (from left to right, so 00112233445566778899aabbccddeeff)
    if (nuc == 'A' || nuc == 'a' || nuc == 0) { // A encoded to 00
        seq[box] &= ~(0b01<<pnt); // encode _0
        seq[box] &= ~(0b10<<pnt); // encode 0_
    }
    else if (nuc == 'C' || nuc == 'c' || nuc == 1) { // C encoded to 01
        seq[box] |= (0b01<<pnt);  // encode _1
        seq[box] &= ~(0b10<<pnt); // encode 0_
    }
    else if (nuc == 'G' || nuc == 'g' || nuc == 2) { // G encoded to 10
        seq[box] &= ~(0b01<<pnt); // encode _0
        seq[box] |= (0b10<<pnt);  // encode 1_
    }
    else if (nuc == 'T' || nuc == 't' || nuc == 3) { // T encoded to 11
        seq[box] |= (0b01<<pnt); // encode _1
        seq[box] |= (0b10<<pnt); // encode 1_
    }
    else { // Randomly encode if not ACGT
        if (std::rand()&1)
            seq[box] &= ~(0b01<<pnt);
        else
            seq[box] |= (0b01<<pnt);
        if (std::rand()&1)
            seq[box] &= ~(0b10<<pnt);
        else
            seq[box] |= (0b10<<pnt);
    }
}

// Read reference sequence and encode into packed form
void read(char* seq_f, std::string& chr, std::uint32_t& len, std::uint32_t* seq) { // (file, chromosomes, length, sequence)
    std::ifstream infile;
    infile.open(seq_f, std::ios::binary);
    infile.seekg(0, std::ios::end);
    std::uint32_t seq_l = infile.tellg();
    char* read(new char[seq_l+1]{});
    infile.seekg(0, std::ios::beg);
    infile.read((char*)read, seq_l);
    infile.close();
    std::uint32_t read_f{0}, nuc_s{0}, nuc_e{0}; // (reading mode, start position, end position)
    std::string chr_name{}; // Chromosome name
    for (std::uint32_t i=0; i<=seq_l; ++i) {
        if (!(read[i])) { // End of file
            if (nuc_e-nuc_s)
                chr += chr_name + '\n' + std::to_string(nuc_e-nuc_s) + '\n'; // Write previous chromosome name and count
            break;
        }
        if (read[i]=='>') { // New entry
            if (nuc_e-nuc_s)
                chr += chr_name + '\n' + std::to_string(nuc_e-nuc_s) + '\n'; // Write previous chromosome name and count
            nuc_s = nuc_e;
            chr_name = "";
            if (read[i+1]=='N' && read[i+2]=='C') // Valid chromosome starting with "NC"
                read_f = 1; // Mode: chromosome name
            else
                read_f = 0; // Mode: skip
        }
        else if (read[i]=='\n' && read_f==1)
            read_f = 2; // Mode: nucleotide
        else if (read_f==1) // Read chromosome name
            chr_name += read[i];
        else if (read_f==2 && read[i]!='\n') { // Read nucleotide
            write_acgt(seq, nuc_e, read[i]);
            nuc_e += 1;
            if ((nuc_e % 1024)==0)
                std::cout << '\r' << "[Read Sequence] " << nuc_e << "/" << seq_l << "                    " << std::flush;
        }
    }
    len = (((nuc_e>>6)+1)<<6); // Padding to 64
    for (std::uint32_t i=nuc_e; i<len; ++i)
        write_acgt(seq, i, 'T');
    delete[] read;
    std::cout << '\r' << "[Read Sequence] " << "Complete" << "                    \n" << std::flush;
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

void build(std::uint32_t len, std::uint32_t* seq, std::uint32_t* sfa, std::uint32_t* bwt, std::uint32_t* occ) {
    const std::uint32_t dig{4};
    const std::uint32_t gid{256};
    std::uint32_t grp[256]{};
    for (std::uint32_t i=0; i<len; ++i)
        grp[nucs(seq, len, i, 4)] += 1;
    std::uint32_t gpc{};
    std::uint32_t cco[4]{};
    for (std::uint32_t g=0; g<256; ++g) {
        std::uint32_t* pmt{new std::uint32_t[grp[g]]{}};
        std::uint32_t* tmp{new std::uint32_t[grp[g]]{}};
        std::uint32_t* pos{new std::uint32_t[grp[g]]{}};
        std::uint32_t* nxt{new std::uint32_t[grp[g]]{}};
        std::uint32_t* cnt{new std::uint32_t[gid]{}};
        std::uint32_t com{0};
        std::uint32_t moc{0};
        std::uint32_t id{};
        for (std::uint32_t i=0; i<len; ++i) {
            if (nucs(seq, len, i, 4)==g) {
                pmt[id] = i;
                pos[id] = 4;
                nxt[id] = grp[g];
                id += 1;
            }
        }
        while (com < grp[g]) {
            if (nxt[com]-com>=2 && nxt[com]-com<=256) {
                std::uint32_t head{pos[com]};
                std::uint32_t num{nxt[com]-com};
                bool flag{1};
                while (flag) {
                    for (std::uint32_t i=1; i<num; ++i) {
                        if (nucs(seq, len, pmt[com]+head, 16)!=nucs(seq, len, pmt[com+i]+head, 16))
                            flag = 0;
                    }
                    if (flag)
                        head += 16;
                }
                for (std::uint32_t i=0; i<num; ++i)
                    pos[com+i] = head;
            }
            for (std::uint32_t i=0; i<gid; ++i)
                cnt[i] = 0;
            for (moc=com; moc<nxt[com]; ++moc) {
                cnt[nucs(seq, len, pmt[moc]+pos[moc], dig)] += 1; 
            }
            for (std::uint32_t j=0; j<cnt[0]; ++j)
                nxt[com+j] = com + cnt[0];
            for (std::uint32_t i=1; i<gid; ++i) {
                cnt[i] += cnt[i-1];
                for (std::uint32_t j=cnt[i-1]; j<cnt[i]; ++j)
                    nxt[com+j] = com + cnt[i];
            }
            for (std::int64_t i=moc-1; i>=com; --i) 
                tmp[com+(--cnt[nucs(seq, len, pmt[i]+pos[i], dig)])] = pmt[i];
            for (std::int64_t i=moc-1; i>=com; --i) {
                pmt[i] = tmp[i];
                pos[i] += dig;
            }
            while (com+1 == nxt[com])
                com += 1;
            std::cout <<'\r' << "[Build Index] " << com+gpc << "/" << len << "                    " << std::flush; 
        }
        for (std::uint32_t i=0; i<grp[g]; ++i) {
            std::uint32_t pre{};
            if (pmt[i]==0) {
                pre = nucs(seq, len, len-1, 1);
                sfa[len/4] = gpc+i;
            }
            else
                pre = nucs(seq, len, pmt[i]-1, 1);
            write_acgt(bwt, gpc+i, (char) pre);
            cco[pre] += 1;
            if ((gpc+i)%4==0)
                sfa[(gpc+i)/4] = pmt[i];
            if ((gpc+i+1)%16==0) {
                for (std::uint32_t j=0; j<4; ++j)
                    occ[((gpc+i)/16)*4+j] = cco[j];
            }
        }
        delete[] pmt;
        delete[] tmp;
        delete[] pos;
        delete[] nxt;
        delete[] cnt;
        gpc += com;
    }
    for (int i=1; i<4; ++i)
        occ[(len>>2)+i] = occ[(len>>2)+i-1] + occ[(((len>>4)-1)<<2)+i-1];
    std::cout << '\r' << "[Build Index] " << "Complete" << "                    \n" << std::flush; 
}

void save(char* seq_f, std::string& chr, std::uint32_t len, std::uint32_t* seq, std::uint32_t* sfa, std::uint32_t* bwt, std::uint32_t* occ) {
    std::cout <<'\r' << "[Save Index] " << "Saving" << "                    " << std::flush; 
    std::ofstream outfile;
    outfile.open(std::string(seq_f)+".chr", std::ios::trunc);
    outfile << chr;
    outfile.close();
    outfile.open(std::string(seq_f)+".seq", std::ios::trunc | std::ios::binary);
    outfile.write((char*) seq, len>>2);
    outfile.close();
    outfile.open(std::string(seq_f)+".sfa", std::ios::trunc | std::ios::binary);
    outfile.write((char*) sfa, len+4);
    outfile.close();
    outfile.open(std::string(seq_f)+".bwt", std::ios::trunc | std::ios::binary);
    outfile.write((char*) bwt, len>>2);
    outfile.close();
    outfile.open(std::string(seq_f)+".occ", std::ios::trunc | std::ios::binary);
    outfile.write((char*) occ, len+16);
    outfile.close();
    std::cout << '\r' << "[Save Index] " << "Complete" << "                    \n" << std::flush; 
}

void index(char** argv) {
    std::time_t start, finish;
    time(&start);
    std::uint32_t len{};
    profile(argv[1], len);
    std::string chr{};
    std::uint32_t* seq{new std::uint32_t[len>>2]{}};
    std::uint32_t* sfa{new std::uint32_t[len+4]{}};
    std::uint32_t* bwt{new std::uint32_t[len>>2]{}};
    std::uint32_t* occ{new std::uint32_t[len+16]{}};
    read(argv[1], chr, len, seq);
    build(len, seq, sfa, bwt, occ);
    save(argv[1], chr, len, seq, sfa, bwt, occ);
    delete[] seq;
    delete[] sfa;
    delete[] bwt;
    delete[] occ;
    time(&finish);
    std::cout << '\r' << "[Finish] Total time: " << difftime(finish, start) << " seconds" << "                    \n" << std::flush;
}

int main(int argc, char** argv) {
    index(argv);
    return 0;
}