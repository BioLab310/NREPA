//
// Created by lin on 21-10-25.
//
#include "operateRef.h"
char* compressRemove(char *mul_ref, uint64_t mul_ref_length,
                     uint64_t *p_remove_len) {
    cout << "开始去除非AGCT字符" << endl;
    //去除非ACGT的字符
    char *p_remove;
    p_remove = (char*) malloc(sizeof(char) * mul_ref_length);
    uint64_t p_remove_length = 0;
    for (uint64_t i = 0; i < mul_ref_length; i++) {
        if (mul_ref[i] == 'A' || mul_ref[i] == 'C' || mul_ref[i] == 'G'
            || mul_ref[i] == 'T') {
            p_remove[p_remove_length++] = mul_ref[i];
        } else {
            continue;
        }
    }
    *p_remove_len = p_remove_length;
    cout << "去除完毕后长度:" << p_remove_length << endl;
    return p_remove;
}

void ReadSeq(char **seq1, uint64_t *seq_length, char *p_ref) {
    uint32_t buffer_size = 256;
    char buffer_line[256];
    //初始化
    memset(buffer_line, 0, buffer_size);

    FILE *fp;
    fp = fopen(p_ref, "r+");
    if (fp == NULL) {
        cout << "file can not be open!" << endl;
        return;
    }

    uint64_t total_size = 0;
    //指针定位到文件尾
    fseek(fp, 0, 2);
    //计算序列总长度
    total_size = ftell(fp);

    char *seq;
    seq = (char*) malloc(sizeof(char) * total_size);

    //指针定位到文件头
    fseek(fp, 0, 0);
    uint64_t len = 0;
    while (fgets(buffer_line, buffer_size - 1, fp) != NULL) {
        if (buffer_line[0] == '>')
            continue;
        else {
            for (uint32_t i = 0; i < buffer_size; i++) {
                if (buffer_line[i] == '\n' || buffer_line[i] == '\0') {
                    break;
                }
                if (buffer_line[i] >= 'a') {
                    buffer_line[i] -= 32;
                }
                seq[len] = buffer_line[i];
                len++;
            }
        }
        memset(buffer_line, 0, buffer_size);
    }
    *seq_length = len;
    *seq1 = seq;
    cout << "the length of seq is: " << len << endl;
}