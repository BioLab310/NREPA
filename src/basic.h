#ifndef _basic_h_
#define _basic_h_
#include <stdint.h>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <cstring>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <algorithm>
#include "stdio.h"
#include <fstream>
#include <ctime>
#include "stdlib.h"
#include <math.h>
#include <iostream>
#include <bitset>
#include <assert.h>
#include <pthread.h>
#include <sys/time.h>
#include <thread>
#include <sstream>
#include <chrono>
#include <unistd.h> 	//access header file
#include <sys/types.h>	//mdkir header file
#include <sys/stat.h> 	//mdkir header file
#include "htslib/sam.h"
using namespace std;
using std::bitset;

#define M 31
struct bit256KmerPara{
    uint32_t kmer1Len;
    uint32_t kmer64Len;
    uint32_t remainer1to64;
    uint64_t codefor1;
};

void kmercpy(uint64_t *dest, const uint64_t * from, uint32_t unitper);
uint32_t cmp256BitKmer(uint64_t*a,uint64_t*b,uint32_t len);
void get_para(struct bit256KmerPara *para1,uint32_t kmer_length);
void ReadSeq(char **seq1,uint32_t *seq_length,const char* p_ref);
char *RefSaveAsSeq(const char *p_ref);
void cal_hash_value_directly_256bit(char *seq,uint64_t * kmer,\
		struct bit256KmerPara para);
void cal_hash_value_indirectly_256bit(char *seq,uint64_t* original,uint64_t* current,\
		struct bit256KmerPara para);
void cal_hash_value_directly_256bit_const(const char *seq,uint64_t * current,\
		struct bit256KmerPara para);
#endif





