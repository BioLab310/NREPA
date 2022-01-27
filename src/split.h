//
// Created by lin on 21-10-25.
//

#ifndef NREPA_SPLIT_H
#define NREPA_SPLIT_H0
#include "basic.h"
#include "compress.h"
#include "load_DBG_full.h"

#endif //NREPA_SPLIT_H
void split_main(int argc, char **argv);
void splitRef(char *p_split_path,uint32_t k_mer);
map<uint64_t,uint64_t> splitNum(char *p_ref_path);
uint64_t res_split(vector<string> &res,map<uint64_t,uint64_t> my_count,char *path);
string split_merge(vector<string> &res,uint32_t k_mer,uint64_t count,uint64_t **start_kmer,uint64_t **end_kmer);