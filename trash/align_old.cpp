#include <algorithm>
#include <cstdint>
#include <ctime>
#include <deque>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

void profile(char* seq_f, std::vector<std::string>& chr_n, std::vector<uint32_t>& chr_c, std::uint32_t& len) {
    std::ifstream infile;
    infile.open(std::string(seq_f)+".chr");
    std::string line;
    while (std::getline(infile, line)) {
        chr_n.push_back(line.substr(0, line.find(' ')));
        std::getline(infile, line);
        chr_c.push_back((std::uint32_t) std::stol(line));
        len += (std::uint32_t) std::stol(line);
    }
    for (int i=1; i<chr_c.size(); ++i)
        chr_c[i] += chr_c[i-1];
    len = (((len>>6)+1)<<6);
    infile.close();
}

void load(char* seq_f, std::uint32_t len, std::uint32_t* seq, std::uint32_t* sfa, std::uint32_t* bwt, std::uint32_t* occ) {
    std::cout <<'\r' << "[Load Index] " << "..." << "                    " << std::flush;
    std::ifstream infile;
    infile.open(std::string(seq_f)+".seq", std::ios::binary);
    infile.read((char*) seq, len>>2);
    infile.close();
    infile.open(std::string(seq_f)+".sfa", std::ios::binary);
    infile.read((char*) sfa, (len>>2)+4);
    infile.close();
    infile.open(std::string(seq_f)+".bwt", std::ios::binary);
    infile.read((char*) bwt, len>>2);
    infile.close();
    infile.open(std::string(seq_f)+".occ", std::ios::binary);
    infile.read((char*) occ, (len>>2)+16);
    infile.close();
    std::cout <<'\r' << "[Load Index] " << "Complete" << "                    " << '\n' << std::flush;
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

std::uint32_t cti(char c) {
    if (c=='A')
        return 0;
    else if (c=='C')
        return 1;
    else if (c=='G')
        return 2;
    else if (c=='T')
        return 3;
    else
        return 0;
}

char itc(std::uint32_t i) {
    if (i==0)
        return 'A';
    else if (i==1)
        return 'C';
    else if (i==2)
        return 'G';
    else if (i==3)
        return 'T';
    else
        return 'A';
}

std::uint32_t lfm(std::uint32_t r, std::uint32_t c, std::uint32_t len, std::uint32_t* sfa, std::uint32_t* bwt, std::uint32_t* occ) {
    std::uint32_t ans{};
    if ((r>>5)&1) {
        ans = occ[(len>>4)+c] + occ[((r>>6)<<2)+c];
        for (std::uint32_t i=r; i<(((r>>6)+1)<<6); ++i) {
            if (nucs(bwt, len, i, 1)==c)
                ans -= 1;
        }
    }
    else {
        if (r<64)
            ans = occ[(len>>4)+c];
        else
            ans = occ[(len>>4)+c] + occ[(((r>>6)-1)<<2)+c];
        for (std::uint32_t i=((r>>6)<<6); i<r; ++i) {
            if (nucs(bwt, len, i, 1)==c)
                ans += 1;
        }
    }
    if (c==3 && r<=sfa[len/16])
        ans += 1;
    return ans;
}

std::uint32_t rpm(std::uint32_t r, std::uint32_t len, std::uint32_t* sfa, std::uint32_t* bwt, std::uint32_t* occ) {
    std::uint32_t row{r};
    std::uint32_t ans{};
    std::uint32_t off{};
    for (off=0; off<len; ++off) {
        if (row%16==0) {
            ans = sfa[row/16];
            break;
        }
        else if (row==sfa[len/16]) {
            ans = 0;
            break;
        }
        row = lfm(row, nucs(bwt, len, row, 1), len, sfa, bwt, occ);
    }
    ans += off;
    return ans;
}

std::string rcseq(std::string& s) {
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

void capitalize(std::string& s) {
    for (std::int32_t i=0; i<s.size(); ++i) {
        if (s[i]=='a')
            s[i] = 'A';
        else if (s[i]=='c')
            s[i] = 'C';
        else if (s[i]=='g')
            s[i] = 'G';
        else if (s[i]=='t')
            s[i] = 'T';
    }
}

void locate(std::string& chr, std::int64_t& pos, std::vector<std::string>& chr_n, std::vector<uint32_t>& chr_c) {
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

struct seed {
    int qs{};
    int qe{};
    std::int64_t rs{};
    std::int64_t re{};
    int fg{};
};

struct cluster {
    int fg{};
    int nc{};
    std::vector<seed> seds;
};

void dp_ed(std::string& sa, std::string& sb, std::vector<int>& aln_i, std::vector<char>& aln_c, int& trc) {
    int dp_i[sa.size()+1][sb.size()+1]{};
    char dp_c[sa.size()+1][sb.size()+1]{};
    std::vector<int> tmp_i;
    std::vector<char> tmp_c;
    for (int i=1; i<=sa.size(); ++i) {
        dp_i[i][0] = i;
        dp_c[i][0] = 'I';
    }
    for (int j=1; j<=sb.size(); ++j) {
        dp_i[0][j] = j;
        dp_c[0][j] = 'D';
    }
    for (int i=1; i<=sa.size(); ++i) {
        for (int j=1; j<=sb.size(); ++j) {
            if (sa[i-1]==sb[j-1]) {
                dp_i[i][j] = dp_i[i-1][j-1];
                dp_c[i][j] = 'M';
                continue;
            }
            else if (dp_i[i-1][j-1]<=dp_i[i][j-1] && dp_i[i-1][j-1]<=dp_i[i-1][j]) {
                dp_i[i][j] = dp_i[i-1][j-1] + 1;
                dp_c[i][j] = 'S';
                continue;
            }
            else if (dp_i[i][j-1]<=dp_i[i-1][j-1] && dp_i[i][j-1]<=dp_i[i-1][j]) {
                dp_i[i][j] = dp_i[i][j-1] + 1;
                dp_c[i][j] = 'D';
                continue;
            }
            else {
                dp_i[i][j] = dp_i[i-1][j] + 1;
                dp_c[i][j] = 'I';
                continue;
            }
        }
    }
    int dx{}, dy{};
    if (trc) {
        trc = 0;
        int tmp{dp_i[sa.size()][0]};
        int j{};
        for (j=1; j<=sb.size(); ++j) {
            if (dp_i[sa.size()][j]<tmp) {
                trc = j;
                tmp = dp_i[sa.size()][j];
            }
        }
        dx = sa.size();
        dy = trc;
    }
    else {
        dx = sa.size();
        dy = sb.size();
    }
    while (dx>0 || dy>0) {
        if (tmp_c.empty() || tmp_c.back()!=dp_c[dx][dy]) {
            tmp_i.push_back(1);
            tmp_c.push_back(dp_c[dx][dy]);
        }
        else
            tmp_i.back() += 1;
        if (dp_c[dx][dy]=='M' || dp_c[dx][dy]=='S') {
            dx -= 1;
            dy -= 1;
        }
        else if (dp_c[dx][dy]=='D')
            dy -= 1;
        else if (dp_c[dx][dy]=='I')
            dx -= 1;
    }
    for (int i=tmp_i.size()-1; i>=0; --i) {
        if (aln_c.empty() || aln_c.back()!=tmp_c[i]) {
            aln_i.push_back(tmp_i[i]);
            aln_c.push_back(tmp_c[i]);
        }
        else
            aln_i.back() += tmp_i[i];
    }
}

void maps(char** argv, std::uint32_t len, std::uint32_t* seq, std::uint32_t* sfa, std::uint32_t* bwt, std::uint32_t* occ, std::vector<std::string>& chr_n, std::vector<std::uint32_t>& chr_c) {
    const int s_wid{10};
    const int gap_1{50};
    const int gap_2{1000};
    const int flk{2};
    std::ifstream qry1_f;
    std::ifstream qry2_f;
    std::ofstream maps_f;
    std::string qry1_n(argv[2]);
    std::string qry2_n(argv[3]);
    std::string maps_n;
    for (int i=0; i<qry1_n.size(); ++i) {
        if (qry1_n[i]!=qry2_n[i])
            break;
        maps_n.push_back(qry1_n[i]);
    }
    maps_n += ((maps_n.back()=='.')?"sam":".sam");
    qry1_f.open(qry1_n);
    qry2_f.open(qry2_n);
    maps_f.open(maps_n, std::ios::trunc);
    for (int i=0; i<chr_n.size(); ++i) {
        maps_f << "@SQ" << '\t';
        maps_f << "SN:" << chr_n[i] << '\t';
        maps_f << "LN:" << (i?(chr_c[i]-chr_c[i-1]):chr_c[0]) << '\n';
    }
    maps_f << "@PG" << '\t' << "ID:vta" << '\t' << "PN:vta" << '\t' << "VN:0.0.1" << '\t';
    maps_f << "CL:" << "./align" << ' ' << argv[1] << ' ' << argv[2] << ' ' << argv[3] << '\n';
    std::vector<std::string> qryn(2);
    std::vector<std::string> qrys(2);
    std::vector<std::string> qals(2);
    std::int64_t qry_c{};
    while (std::getline(qry1_f, qryn[0])) {
        std::getline(qry1_f, qrys[0]);
        std::getline(qry1_f, qals[0]);
        std::getline(qry1_f, qals[0]);
        std::getline(qry2_f, qryn[1]);
        std::getline(qry2_f, qrys[1]);
        std::getline(qry2_f, qals[1]);
        std::getline(qry2_f, qals[1]);
        qryn[0] = qryn[0].substr(1, qryn[0].find(' ')-1);
        qryn[1] = qryn[1].substr(1, qryn[1].find(' ')-1);
        std::vector<std::string> seqs{qrys[0], rcseq(qrys[0]), qrys[1], rcseq(qrys[1])};
        std::vector<seed> seds;
        std::deque<std::int64_t> sedt(s_wid);
        for (int fg=0; fg<4; ++fg) {
            std::string qry = seqs[fg];
            for (int qe=seqs[fg].size(); qe>0; --qe) {
                std::int64_t head = 0;
                std::int64_t tail = len;
                for (int qp=qe-1; qp>=0; --qp) {
                    head = lfm(head, cti(qry[qp]), len, sfa, bwt, occ);
                    tail = lfm(tail, cti(qry[qp]), len, sfa, bwt, occ);
                    if ((tail-head)<=s_wid) {
                        for (std::int64_t rps=head; rps<tail; ++rps) {
                            std::int64_t rp = rpm(rps, len, sfa, bwt, occ);
                            if (std::find(sedt.begin(), sedt.end(), rp+qry.size()-qp)==sedt.end()) {
                                int qs{}; 
                                for (qs=qp-1; qs>=0; --qs) {
                                    if (qry[qs]!=itc(nucs(seq, len, rp-qp+qs, 1)))
                                        break;
                                }
                                qs += 1;
                                seds.push_back({qs, qe, rp-qp+qs, rp+qe-qp, fg});
                            }
                            sedt.push_back(rp+qry.size()-qp);
                        }
                        for (int i=0; i<s_wid-(tail-head); ++i)
                            sedt.push_back(0);
                        for (int i=0; i<s_wid; ++i)
                            sedt.pop_front();
                        break;
                    }
                }
            }
        }
        std::sort(seds.begin(), seds.end(), [](seed a, seed b){return (a.rs<b.rs) ? 1 : 0;});
        std::vector<cluster> cls;
        int now[4]{-1, -1, -1, -1};
        int chs[4]{-1, -1, -1, -1};
        int val[4]{};
        int l{}, fg{};
        for (int i=0; i<seds.size(); ++i) {
            fg = seds[i].fg;
            if ((now[fg])==-1) {
                cls.push_back({seds[i].fg, seds[i].qe-seds[i].qs, {{seds[i].qs, seds[i].qe, seds[i].rs, seds[i].re, seds[i].fg}}});
                now[fg] = cls.size()-1;
                continue;
            }
            if ((seds[i].rs > cls[now[fg]].seds.back().re + gap_1)) {
                if (cls[now[fg]].nc > val[fg]) {
                    chs[fg] = now[fg];
                    val[fg] = cls[now[fg]].nc;
                }
                for (int j=now[fg]-1; j>=0; --j) {
                    if (cls[now[fg]].seds.back().rs>cls[j].seds.back().rs+gap_2)
                        break;
                    if (fg>>1 != cls[j].fg>>1) {
                        if (cls[now[fg]].nc+cls[j].nc>val[fg]) {
                            chs[fg] = now[fg];
                            val[fg] = cls[now[fg]].nc+cls[j].nc;
                        }
                        if (cls[now[fg]].nc+cls[j].nc>val[cls[j].fg]) {
                            chs[cls[j].fg] = j;
                            val[cls[j].fg] = cls[now[fg]].nc+cls[j].nc;
                        }
                    }
                }
                cls.push_back({seds[i].fg, seds[i].qe-seds[i].qs, {{seds[i].qs, seds[i].qe, seds[i].rs, seds[i].re, seds[i].fg}}});
                now[fg] = cls.size()-1;
            }
            else {
                if (seds[i].rs<cls[now[fg]].seds.back().re) {
                    if ((seds[i].qe-seds[i].qs)>(cls[now[fg]].seds.back().qe-cls[now[fg]].seds.back().qs))
                        cls[now[fg]].seds.pop_back();
                    else
                        continue;
                }
                else if (seds[i].qs<cls[now[fg]].seds.back().qe) {
                    if ((seds[i].qe-seds[i].qs)>(cls[now[fg]].seds.back().qe-cls[now[fg]].seds.back().qs))
                        cls[now[fg]].seds.pop_back();
                    else
                        continue;
                }
                if (cls[now[fg]].seds.empty()) {
                    cls[now[fg]].seds.push_back({seds[i].qs, seds[i].qe, seds[i].rs, seds[i].re, seds[i].fg});
                    cls[now[fg]].nc = seds[i].qe-seds[i].qs;
                }
                else {
                    l = std::min(std::min((std::int64_t) seds[i].qe-cls[now[fg]].seds.back().qe, seds[i].re-cls[now[fg]].seds.back().re), (std::int64_t) seds[i].qe-seds[i].qs);
                    if (l>0) {
                        cls[now[fg]].nc += l;
                        cls[now[fg]].seds.push_back({seds[i].qe-l, seds[i].qe, seds[i].re-l, seds[i].re, seds[i].fg});
                    }                
                }
            }
        }
        for (int i=0; i<4; ++i) {
            if (now[i]==-1)
                continue;
            if (cls[now[i]].nc > val[i]) {
                chs[i] = now[i];
                val[i] = cls[now[i]].nc;
            }
            for (int j=now[i]-1; j>=0; --j) {
                if (cls[now[i]].seds.back().rs>cls[j].seds.back().rs+gap_2)
                    break;
                if (i>>1 != cls[j].fg>>1) {
                    if (cls[now[i]].nc+cls[j].nc>val[i]) {
                        chs[i] = now[i];
                        val[i] = cls[now[i]].nc+cls[j].nc;
                    }
                    if (cls[now[i]].nc+cls[j].nc>val[cls[j].fg]) {
                        chs[cls[j].fg] = j;
                        val[cls[j].fg] = cls[now[i]].nc+cls[j].nc;
                    }
                }
            }
        }
        std::string chr[2]{};
        std::int64_t pos[2]{};
        int rf[2]{};
        std::vector<int> aln_i[2];
        std::vector<char> aln_c[2];
        for (int q=0; q<2; ++q) {
            if (chs[q*2+1]==-1 && chs[q*2]==-1) {
                chr[q] = "NA";
                pos[q] = 0;
                continue;
            }
            if (chs[q*2+1]==-1 || chs[q*2]==-1)
                rf[q] = (chs[q*2+1]==-1)?0:1;
            else if (val[q*2+1]>val[q*2])
                rf[q] = 1;
            else if (val[q*2+1]==val[q*2] && cls[chs[q*2+1]].nc>cls[chs[q*2]].nc)
                rf[q] = 1;
            else
                rf[q] = 0;
            std::string qry = seqs[q*2+rf[q]];
            std::vector<seed> sed_m = cls[chs[q*2+rf[q]]].seds;
            std::string sa;
            std::string sb;
            sa = qry.substr(0, sed_m.front().qs);
            std::reverse(sa.begin(), sa.end());
            sb = "";
            for (int i=0; i<sa.size()*flk; ++i)
                sb.push_back(itc(nucs(seq, len, sed_m.front().rs-1-i, 1)));
            int trc{1};
            dp_ed(sa, sb, aln_i[q], aln_c[q], trc);
            pos[q] = sed_m.front().rs - trc;
            std::reverse(aln_i[q].begin(), aln_i[q].end());
            std::reverse(aln_c[q].begin(), aln_c[q].end());
            if (aln_c[q].empty() || aln_c[q].back()!='M') {
                aln_i[q].push_back(sed_m.front().qe-sed_m.front().qs);
                aln_c[q].push_back('M');
            }
            else
                aln_i[q].back() += sed_m.front().qe-sed_m.front().qs;
            for (int s=1; s<sed_m.size(); ++s) {
                sa = qry.substr(sed_m[s-1].qe, sed_m[s].qs-sed_m[s-1].qe);
                sb = "";
                for (std::int64_t i=sed_m[s-1].re; i<sed_m[s].rs; ++i)
                    sb.push_back(itc(nucs(seq, len, i, 1)));
                trc = 0;
                dp_ed(sa, sb, aln_i[q], aln_c[q], trc);
                if (aln_c[q].empty() || aln_c[q].back()!='M') {
                    aln_i[q].push_back(sed_m[s].qe-sed_m[s].qs);
                    aln_c[q].push_back('M');
                }
                else
                    aln_i[q].back() += sed_m[s].qe-sed_m[s].qs;
            }
            sa = qry.substr(sed_m.back().qe, qry.size()-sed_m.back().qe);
            sb = "";
            for (int i=0; i<sa.size()*flk; ++i)
                sb.push_back(itc(nucs(seq, len, sed_m.back().re+i, 1)));
            trc = 1;
            dp_ed(sa, sb, aln_i[q], aln_c[q], trc);
            locate(chr[q], pos[q], chr_n, chr_c);
        }
        for (int q=0; q<2; ++q) {
            maps_f << qryn[q] << '\t';
            int flag{1};
            if (chr[q]==chr[1-q] && pos[q]-pos[1-q]<gap_2 && pos[1-q]-pos[q]<gap_2 && rf[q]!=rf[1-q])
                flag += 2;
            if (pos[q]==0)
                flag += 4;
            if (pos[1-q]==0)
                flag += 8;
            if (rf[q]==1)
                flag += 16;
            if (rf[1-q]==1)
                flag += 32;
            if (q==0)
                flag += 64;
            if (q==1)
                flag += 128;
            maps_f << flag << '\t';
            maps_f << chr[q] << '\t';
            maps_f << pos[q] << '\t';
            maps_f << 255 << '\t';
            std::string cigar;
            std::vector<int> nla_i;
            std::vector<char> nla_c;
            for (int i=0; i<aln_i[q].size(); ++i) {
                if (aln_c[q][i]=='M' || aln_c[q][i]=='S') {
                    if (i==0 || (aln_c[q][i-1]!='M' && aln_c[q][i-1]!='S')) {
                        nla_i.push_back(aln_i[q][i]);
                        nla_c.push_back('M');
                    }
                    else {
                        nla_i.back() += aln_i[q][i];
                    }
                }
                else {
                    nla_i.push_back(aln_i[q][i]);
                    nla_c.push_back(aln_c[q][i]);
                }
            }
            if (pos[q]==0)
                maps_f << "NA";
            for (int i=0; i<nla_i.size(); ++i)
                maps_f << nla_i[i] << nla_c[i];
            maps_f << '\t';
            if (chr[1-q]==chr[q])
                maps_f << "=" << '\t';
            else
                maps_f << chr[1-q] << '\t';
            maps_f << pos[1-q] << '\t';
            int64_t tmpl{};
            if (flag & 2) {
                int qt = (pos[1]>pos[0])?1:0;
                tmpl = pos[qt] - pos[1-qt];
                for (int i=0; i<aln_i[qt].size(); ++i) {
                    if (aln_c[qt][i]!='I')
                        tmpl += aln_i[qt][i];
                }
                if (rf[q])
                    tmpl *= -1;
            }
            maps_f << tmpl << '\t';
            maps_f << qrys[q] << '\t';
            maps_f << qals[q] << '\t';
            maps_f << "aln:";
            for (int i=0; i<aln_i[q].size(); ++i)
                maps_f << aln_c[q][i] << aln_i[q][i];
            maps_f << '\n';
        }
        std::cout <<'\r' << "[Align Read] " << ++qry_c << "                    " << std::flush;
    }
    qry1_f.close();
    qry2_f.close();
    maps_f.close();
    std::cout << '\r' << "[Align Reads] " << "Complete" << "                    \n" << std::flush; 
}

void align(char** argv) {
    std::time_t start, finish;
    time(&start);
    std::vector<std::string> chr_n;
    std::vector<std::uint32_t> chr_c;
    std::uint32_t len{};
    profile(argv[1], chr_n, chr_c, len);
    std::uint32_t* seq{new std::uint32_t[len>>4]{}};
    std::uint32_t* sfa{new std::uint32_t[(len>>4)+1]{}};
    std::uint32_t* bwt{new std::uint32_t[len>>4]{}};
    std::uint32_t* occ{new std::uint32_t[(len>>4)+4]{}};
    load(argv[1], len, seq, sfa, bwt, occ);
    maps(argv, len, seq, sfa, bwt, occ, chr_n, chr_c);
    delete[] seq;
    delete[] sfa;
    delete[] bwt;
    delete[] occ;
    time(&finish);
    std::cout << '\r' << "[Finish] Total time: " << difftime(finish, start) << " seconds" << "                    \n" << std::flush;
}

int main(int argc, char** argv) {
    align(argv);
    return 0;
}