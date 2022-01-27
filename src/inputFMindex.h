/*
 * inputFMindex.h
 *
 *  Created on: Apr 15, 2021
 *      Author: bio
 */

#ifndef INPUTFMINDEX_H_
#define INPUTFMINDEX_H_

#include "basic.h"
struct sFMindex
{
	uint32_t sa_gap;
	uint32_t occ_gap;
	char *b;
	uint32_t *c;
	uint32_t *occa;
	uint32_t *sa;
	uint32_t occ_num;
	uint32_t sa_num;
};

extern struct sFMindex gnFMidx,grFMidx;
void *read_binfile(char *path);
void read_bfile2index(const char *path, sFMindex &FMidx, uint32_t num);
void free_FMindex(sFMindex &FMidx);

#endif /* INPUTFMINDEX_H_ */
