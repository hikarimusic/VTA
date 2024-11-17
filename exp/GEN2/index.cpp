#include <iostream>
#include <fstream>
#include <cstring>
#include <string>

void read(char* filename, uint32_t& len, std::string& chr, uint8_t* str) {
	// Read whole file
	std::ifstream file(filename, std::ios::binary);
	file.seekg(0, std::ios::end);
	len = file.tellg();
	file.seekg(0, std::ios::beg);
	file.read(reinterpret_cast<char*>(str), len);
	file.close();
    
	// In place processing
	uint32_t seq_start = 0;
	uint32_t read_pos = 0;
	uint32_t write_pos = 0;
	std::string header = "";
	for (; read_pos<len; ++read_pos) {
    	char c = str[read_pos];
    	if (c == '>') {
        	if (header.size()>2 && header.substr(0, 2)=="NC") {
            	header = header.substr(0, header.find(" "));
            	chr += (header + " " + std::to_string(write_pos - seq_start) + "\n");
            	seq_start = write_pos;
        	}
        	header.clear();
        	read_pos++;
        	while (read_pos < len && str[read_pos] != '\n') {
            	header += str[read_pos];
            	read_pos++;
        	}
        	continue;
    	}
    	if (c == 'A' || c == 'a')
        	str[write_pos++] = 'A';
    	else if (c == 'C' || c == 'c')
        	str[write_pos++] = 'C';
    	else if (c == 'G' || c == 'g')
        	str[write_pos++] = 'G';
    	else if (c == 'T' || c == 't')
        	str[write_pos++] = 'T';
    	else if (c != '\n') {
        	str[write_pos++] = 'Z';
    	}
    	if ((read_pos % 1000000)==0)
        	std::cout << '\r' << "[Read Sequence] " << read_pos << "/" << len << "                	" << std::flush;
	}
	if (header.size()>2 && header.substr(0, 2)=="NC") {
    	header = header.substr(0, header.find(" "));
    	chr += (header + " " + std::to_string(write_pos - seq_start) + "\n");
	}
	str[write_pos++] = '$';
	len = write_pos;
	std::cout << '\r' << "[Read Sequence] " << "Complete" << "                	\n" << std::flush;
}

const uint8_t mask[]={0x80, 0x40, 0x20, 0x10, 0x08, 0x04, 0x02, 0x01};
#define min1 4294967295
#define tget(i) ((T[(i)/8] & mask[(i)%8]) ? 1 : 0)
#define tset(i, b) T[(i)/8] = (b) ? (mask[(i)%8] | T[(i)/8]) : ((~mask[(i)%8]) & T[(i)/8])
#define chr(i) (cs==sizeof(uint32_t) ? ((uint32_t*)S)[i] : ((uint8_t*)S)[i])
#define isLMS(i) (i>0 && tget(i) && !tget(i-1))

// find the start or end of each bucket
inline void getBuckets(uint8_t* S, uint32_t* B, uint32_t n, uint32_t k, uint32_t cs, bool end) {
	uint32_t sum=0;
	// clear all buckets
	for (uint32_t i=0; i<=k; ++i)
    	B[i]=0;
	// compute the size of each bucket
	for (uint32_t i=0; i<n; ++i)
    	++B[chr(i)];
	for (uint32_t i=0; i<=k; ++i) {
    	sum += B[i];
    	B[i] = end ? sum : sum-B[i];
	}
}

// compute SAl
inline void induceSAl(uint8_t* T, uint32_t* SA, uint8_t* S, uint32_t* B, uint32_t n, uint32_t k, uint8_t cs, bool end) {
	getBuckets(S, B, n, k, cs, end);
	for (uint32_t i=0; i<n; ++i) {
    	if (SA[i]==min1 || SA[i]==0)
        	continue;
    	uint32_t j = SA[i]-1;
    	if (!tget(j))
        	SA[B[chr(j)]++] = j;
	}
}
// compute SAs
inline void induceSAs(uint8_t* T, uint32_t* SA, uint8_t* S, uint32_t* B, uint32_t n, uint32_t k, uint8_t cs, bool end) {
	getBuckets(S, B, n, k, cs, end);
	for (uint32_t i1=n; i1>=1; --i1) {
    	uint32_t i = i1 - 1;
    	if (SA[i]==min1 || SA[i]==0)
        	continue;
    	uint32_t j = SA[i]-1;
    	if (tget(j))
        	SA[--B[chr(j)]] = j;
	}
}

void SA_IS(uint8_t* S, uint32_t* SA, uint32_t* B, uint8_t* T, uint32_t n, uint32_t k, int8_t cs) {
	// LS-type array in bits
	tset(n-1, 1);
	tset(n-2, 0);
	for (uint32_t i1=n-2; i1>=1; --i1) {
    	uint32_t i = i1 - 1;
    	tset(i, (chr(i)<chr(i+1) || (chr(i)==chr(i+1) && tget(i+1)==1)) ? 1 : 0);
	}
	// stage 1: reduce the problem by at least 1/2
	// sort all the S-substrings
	getBuckets(S, B, n, k, cs, true);
	for (uint32_t i=0; i<n; ++i)
    	SA[i] = min1;
	for (uint32_t i=1; i<n; ++i) {
    	if (isLMS(i))
        	SA[--B[chr(i)]] = i;
	}
	induceSAl(T, SA, S, B, n, k, cs, false);
	induceSAs(T, SA, S, B, n, k, cs, true);
	// compact all the sorted substrings into the first n1 items of SA
	uint32_t n1=0;
	for (uint32_t i=0; i<n; ++i) {
    	if (isLMS(SA[i]))
        	SA[n1++]=SA[i];
	}
	// find the lexicographic names of substrings
	for (uint32_t i=n1; i<n; ++i)
    	SA[i] = min1;
	uint32_t name=0, prev=min1;
	for (uint32_t i=0; i<n1; ++i) {
    	uint32_t pos = SA[i];
    	bool diff = false;
    	for (uint32_t d=0; d<n; ++d) {
        	if (prev==min1 || chr(pos+d)!=chr(prev+d) || tget(pos+d)!=tget(prev+d)) {
            	diff = true;
            	break;
        	}
        	else if (d>0 && (isLMS(pos+d) || isLMS(prev+d)))
            	break;
    	}
    	if (diff) {
        	++name;
        	prev = pos;
    	}
    	pos = (pos%2==0) ? pos/2 : (pos-1)/2;
    	SA[n1+pos] = name-1;
	}
	for (uint32_t i=n-1, j=n-1; i>=n1; --i) {
    	if (SA[i]!=min1)
        	SA[j--] = SA[i];
	}
	// stage 2: solve the reduced problem
	uint32_t *SA1=SA, *s1=SA+n-n1, *B1=B;
	uint8_t *T1 = T+n/8+1;
	if (name<n1) {
    	SA_IS((uint8_t*)s1, SA1, B1, T1, n1, name-1, sizeof(uint32_t));
	}
	else {
    	for (uint32_t i=0; i<n1; ++i)
        	SA1[s1[i]] = i;
	}
	// stage 3: induce the result for the original problem
	// put all the LMS characters into their buckets
	getBuckets(S, B, n, k, cs, true);
	for (uint32_t i=1, j=0; i<n; ++i) {
    	if (isLMS(i))
        	s1[j++] = i;
	}
	for (uint32_t i=0; i<n1; ++i)
    	SA1[i] = s1[SA1[i]];
	for (uint32_t i=n1; i<n; ++i)
    	SA[i] = min1;
	for (uint32_t i1=n1; i1>=1; --i1) {
    	uint32_t i = i1 - 1;
    	uint32_t j = SA[i];
    	SA[i] = min1;
    	SA[--B[chr(j)]] = j;
	}
	induceSAl(T, SA, S, B, n, k, cs, false);
	induceSAs(T, SA, S, B, n, k, cs, true);
}

void index(char** argv) {
	uint32_t len = 0;
	std::string chr;
	uint8_t* str = new uint8_t[4000000000]; // 4G
	uint32_t* sfa = new uint32_t[4000000000]; // 16G
	uint32_t* bkt = new uint32_t[2000000000]; // 8G
	uint8_t* typ = new uint8_t[1000000000]; // 1G
    
	// Read sequence
	read(argv[1], len, chr, str);
    
	// Construct suffix array
	SA_IS(str, sfa, bkt, typ, len, 255, 1);
    
	delete[] str;
	delete[] sfa;
	delete[] bkt;
	delete[] typ;
}

int main(int argc, char** argv) {
	index(argv);
	return 0;
}