#include <cstdlib>
#include <cstdint>
#include <ctime>
#include <fstream>
#include <iostream>
#include <string>

void write_acgt(std::uint32_t* seq, std::uint32_t pos, char nuc) {
    std::uint32_t box = pos >> 4;
    std::uint32_t pnt = (0b1111 - pos & 0b1111) << 1;
    if (nuc == 'A' || nuc == 'a' || nuc == 0) {
        seq[box] &= ~(0b01<<pnt);
        seq[box] &= ~(0b10<<pnt);
    }
    else if (nuc == 'C' || nuc == 'c' || nuc == 1) {
        seq[box] |= (0b01<<pnt);
        seq[box] &= ~(0b10<<pnt);
    }
    else if (nuc == 'G' || nuc == 'g' || nuc == 2) {
        seq[box] &= ~(0b01<<pnt);
        seq[box] |= (0b10<<pnt);
    }
    else if (nuc == 'T' || nuc == 't' || nuc == 3) {
        seq[box] |= (0b01<<pnt);
        seq[box] |= (0b10<<pnt);
    }
    else {
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

void read(char* seq_f, std::string& chr, std::uint32_t& len, std::uint32_t* seq) {
    std::ifstream infile;
    infile.open(seq_f, std::ios::binary);
    infile.seekg(0, std::ios::end);
    std::uint32_t seq_l = infile.tellg();
    char* read(new char[seq_l+1]{});
    infile.seekg(0, std::ios::beg);
    infile.read((char*)read, seq_l);
    infile.close();
    std::uint32_t read_f{0}, nuc_s{0}, nuc_e{0};
    std::string chr_name{};
    for (std::uint32_t i=0; i<=seq_l; ++i) {
        if (!(read[i])) {
            if (nuc_e-nuc_s)
                chr += chr_name + '\n' + std::to_string(nuc_e-nuc_s) + '\n';
            break;
        }
        if (read[i]=='>') {
            if (nuc_e-nuc_s) 
                chr += chr_name + '\n' + std::to_string(nuc_e-nuc_s) + '\n';
            nuc_s = nuc_e;
            chr_name = "";
            if (read[i+1]=='N' && read[i+2]=='C')
                read_f = 1;
            else
                read_f = 0;
        }
        else if (read[i]=='\n' && read_f==1)
            read_f = 2;
        else if (read_f==1)
            chr_name += read[i];
        else if (read_f==2 && read[i]!='\n') {
            write_acgt(seq, nuc_e, read[i]);
            nuc_e += 1;
            if ((nuc_e % 1024)==0)
                std::cout << '\r' << "[Read Sequence] " << nuc_e << "/" << seq_l << "                    " << std::flush;
        }
    }
    len = (((nuc_e>>6)+1)<<6);
    for (std::uint32_t i=nuc_e; i<len; ++i)
        write_acgt(seq, i, 'T');
    delete[] read;
    std::cout << '\r' << "[Read Sequence] " << "Complete" << "                    \n" << std::flush;
} 

std::uint32_t nucs(std::uint32_t* seq, std::uint32_t len, std::uint32_t p, std::uint32_t r) {
    std::uint32_t b = p >> 4;
    std::uint32_t h = p & 0b1111;
    std::uint32_t t = h + r;
    if (p >= len) {
        return 0;
    }
    else if (p+r >= len) {
        t -= 16;
        return (((seq[b]<<(h<<1))>>(h<<1))<<(t<<1));
    }
    else if (t <= 16) {
        return ((seq[b]<<(h<<1))>>((16-r)<<1));
    }
    else {
        t -= 16;
        return  ((((seq[b]<<(h<<1))>>(h<<1))<<(t<<1))+(seq[b+1]>>((16-t)<<1)));
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
            if (com % 1024 == 0)
                std::cout <<'\r' << "[Build Index] " << com+gpc << "/" << len << "                    " << std::flush; 
        }
        for (std::uint32_t i=0; i<grp[g]; ++i) {
            std::uint32_t pre{};
            if (pmt[i]==0) {
                pre = nucs(seq, len, len-1, 1);
                sfa[len/16] = gpc+i;
            }
            else
                pre = nucs(seq, len, pmt[i]-1, 1);
            write_acgt(bwt, gpc+i, (char) pre);
            cco[pre] += 1;
            if ((gpc+i)%16==0)
                sfa[(gpc+i)/16] = pmt[i];
            if ((gpc+i+1)%64==0) {
                for (std::uint32_t j=0; j<4; ++j)
                    occ[((gpc+i)/64)*4+j] = cco[j];
            }
        }
        delete[] pmt;
        delete[] tmp;
        delete[] pos;
        delete[] cnt;
        gpc += com;
    }
    for (int i=1; i<4; ++i)
        occ[(len>>4)+i] = occ[(len>>4)+i-1] + occ[(((len>>6)-1)<<2)+i-1];
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
    outfile.write((char*) sfa, (len>>2)+4);
    outfile.close();
    outfile.open(std::string(seq_f)+".bwt", std::ios::trunc | std::ios::binary);
    outfile.write((char*) bwt, len>>2);
    outfile.close();
    outfile.open(std::string(seq_f)+".occ", std::ios::trunc | std::ios::binary);
    outfile.write((char*) occ, (len>>2)+16);
    outfile.close();
    std::cout << '\r' << "[Save Index] " << "Complete" << "                    \n" << std::flush; 
}

void index(char** argv) {
    std::time_t start, finish;
    time(&start);
    std::string chr{};
    std::uint32_t len{};
    std::uint32_t* seq{new std::uint32_t[268435456]{}};
    read(argv[1], chr, len, seq);
    std::uint32_t* sfa{new std::uint32_t[268435456]{}};
    std::uint32_t* bwt{new std::uint32_t[268435456]{}};
    std::uint32_t* occ{new std::uint32_t[268435456]{}};
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