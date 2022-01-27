//
// Created by lin on 21-10-25.
//

#ifndef NREPA_UTILS_H
#define NREPA_UTILS_H

#include "basic.h"
#include "operateRef.h"
#include "load_DBG_full.h"
#include "split.h"
#include "compress.h"
#define MAXSIZE 10000000000
void utils_main(int argc, char **argv);

void split_detail(int argc, char **argv);
void merge(uint32_t f_count,char **p_ref_path);
void copy_utils(uint64_t copy_size);
void split_compress_ref(char *p_ref_path,uint64_t start,uint64_t span,char *p_sav_path);
void splitNum_v2(char *p_ref_path);
void readArray(char *outputFileName);
void writeArray();

#endif //NREPA_UTILS_H
