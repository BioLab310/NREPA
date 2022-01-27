//
// Created by lin on 21-10-25.
//

#ifndef NREPA_NRCPA_H
#define NREPA_NRCPA_H

#include "basic.h"
#include "operateRef.h"
#include "Hash.h"
#include "Binary_Search.h"
#include "load_DBG_full.h"

//保存不同长度的k_mer,k_code_size也不同
extern uint32_t k_code_size;

struct para_multi{
    char **p_ref;
    vector<bool> *is_find;
    uint64_t start;
    uint64_t end;
    uint32_t k_code_size;
    uint32_t k_mer;
    uint64_t **read;
    uint32_t ***arrptr;
};

void compress_main(int argc, char **argv);

void compress_by_array(char *mul_ref,uint64_t mul_ref_length,uint32_t k_mer,uint64_t start,uint64_t end,\
    char *new_sort_file);

void toArray(char *mul_ref,uint64_t mul_ref_length,uint32_t k_mer,uint64_t start,uint64_t end,uint32_t thread_num,\
    char *sort_file_path);

void *Parall_compress_ref(void *p);

void saveToFile(char *mul_ref,uint64_t mul_ref_length,uint32_t k_mer,uint64_t start,uint64_t end,\
    char *new_sort_file);

//决定方法的入口
void compressing(char *mul_ref, uint64_t mul_ref_length, uint32_t k_mer, uint64_t start, uint64_t end, uint32_t thread_num, \
    char *sort_file_path, char *new_sort_file, uint32_t fun_opt);

uint32_t kmer_codeSize(uint32_t k_mer);

#endif //NREPA_NRCPA_H
