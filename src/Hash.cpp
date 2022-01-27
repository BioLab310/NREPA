/*
 * Hash.cpp
 *
 *  Created on: Jan 24, 2019
 *      Author: bio
 */
#include "Hash.h"

uint32_t bit256hashFFunction(uint64_t *a,struct bit256KmerPara para)
{
	uint64_t tmp=0;
	if(para.kmer64Len>1)
	{
		tmp=tmp|(a[para.kmer64Len-2]<<para.remainer1to64);
		tmp=tmp|a[para.kmer64Len-1];
	}
	else
	{
		tmp=a[para.kmer64Len-1];
	}
	return tmp&(0xFFFFFF);
}
struct NodeBit** bit256initialHashFTable()
{
	struct NodeBit** pk;
	pk=(struct NodeBit**)malloc(sizeof(struct NodeBit*)*HashFSize);
	for(uint32_t i=0;i<HashFSize;i++)
	{
		pk[i]=NULL;
	}
	return pk;
}
struct NodeBit * bit256initialSingleHashFTable()
{
	struct NodeBit* p;
	p=(struct NodeBit*)malloc(sizeof(struct NodeBit));
	Node_initial_bit(p);
	return p;
}
void bit256insertHashFTable(struct NodeBit** pk,struct nodeBit c_tmp_hashtable,struct bit256KmerPara para)
{
	uint32_t tmp=bit256hashFFunction(c_tmp_hashtable.hashValue,para);
	if(pk[tmp]==NULL)
	{
		pk[tmp]=bit256initialSingleHashFTable();
	}
	Insert_Value_bit(pk+tmp,c_tmp_hashtable,para);
}
int64_t getHashFTableValue(struct NodeBit** pk,uint64_t *a,struct bit256KmerPara para)
{
	uint32_t tmp=bit256hashFFunction(a,para);
	if(pk[tmp]!=NULL)
	{
		return MappingHashValueToID_bit(pk[tmp],a,para);
	}
	return -1;
}
void bit256freeHashFTable(struct NodeBit **ph,struct bit256KmerPara para)
{
	for(uint32_t i=0;i<HashFSize;i++)
	{
		if(ph[i]!=NULL)
		{
			destory_tree_bit(ph[i],para);
		}
	}
	free(ph);
}


