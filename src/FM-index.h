//
// Created by lin on 2021/11/18.
//

#ifndef NREPA_FM_INDEX_H
#define NREPA_FM_INDEX_H
#include "basic.h"
#include "B+tree.h"
struct build_para
{
    uint32_t level;
    char * ref_path;
    uint32_t thread_num;
    uint32_t sa_gap;
    uint32_t occ_gap;
    uint64_t max_len;
};
void fm_main(int argc, char **argv);
void build(struct build_para &para_build);
void build_occA(struct build_para &para_build);
#endif //NREPA_FM_INDEX_H
