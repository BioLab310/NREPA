/*
 * Hash.h
 *
 *  Created on: Jan 24, 2019
 *      Author: bio
 */

#ifndef HASH_H_
#define HASH_H_

#include "basic.h"
#include "BplusTreeBit.h"
#define HashFSize 0x1000000
uint32_t bit256hashFFunction(uint64_t *a,\
		struct bit256KmerPara para);
struct NodeBit** bit256initialHashFTable();
struct NodeBit * bit256initialSingleHashFTable();
void bit256insertHashFTable(struct NodeBit** pk,\
		struct nodeBit c_tmp_hashtable,\
		struct bit256KmerPara para);
void bit256freeHashFTable(struct NodeBit **ph,\
		struct bit256KmerPara para);
int64_t getHashFTableValue(struct NodeBit** pk,\
		uint64_t *a,\
		struct bit256KmerPara para);
#endif /* HASH_H_ */
