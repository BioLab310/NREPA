//
// Created by lin on 21-10-25.
//

#ifndef NREPA_OPERATEREF_H
#define NREPA_OPERATEREF_H

#include "basic.h"
//去除非ACGT的字符
char *compressRemove(char *mul_ref, uint64_t mul_ref_length,
                     uint64_t *p_remove_len);

//读基因组文件
void ReadSeq(char **seq1, uint64_t *seq_length, char *p_ref);

#endif //NREPA_OPERATEREF_H
