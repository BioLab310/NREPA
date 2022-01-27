/*
 * BPlusTree_full.h
 *
 *  Created on: Sep 25, 2019
 *      Author: pluto
 */

#ifndef BINARY_SEARCH_H_
#define BINARY_SEARCH_H_

#include "basic.h"
#include <limits.h>
#define B 31

uint32_t  **generate_array(uint32_t total_num, uint32_t save2file);
void free_genarray(uint32_t  ***arrptr);
uint32_t find_arrindex(uint32_t  **arrptr, uint64_t* const kmer_arr,uint64_t kmer_inf);
uint32_t find_arrindexN(uint32_t  **arrptr, uint64_t* const kmer_arr,uint64_t *kmer_inf,uint32_t len);
uint32_t find_supathindex(uint32_t  **arrptr, uint32_t *super_arr,uint32_t pos);


template <typename T>
T  **Tgenerate_array(T total_num)
{
	if(total_num <= 0 || total_num < B)
	{
		printf("Wrong Input Number : Tgenerate_array!\n");
		exit(-1);
	}
	T level_num;
	float rslt;
	rslt = log10(total_num) / log10(B);
	level_num = (T)ceil(rslt);
	T ** arrptr;
	arrptr = new T*[level_num];
	arrptr[0] = new T;
	arrptr[0][0] = level_num;
	T remain_num,current_num;
	current_num = total_num;
	for(T i = level_num; i > 1; i--)
	{
		if(i == level_num)
		{
			remain_num = current_num % B;
		}
		current_num = (T)ceil((double)current_num / B);
		arrptr[i-1] = new T[current_num+1];
		arrptr[i-1][0] = current_num;
		for(T j = 1; j <= current_num; j++)
		{
			if(i == level_num)
			{
				if(j == current_num)
				{
					arrptr[i-1][j] = total_num;
					break;
				}
				arrptr[i-1][j] = j * B;
			}
			else
			{
				if(j == current_num)
				{
					arrptr[i-1][j] = total_num;
					break;
				}
				arrptr[i-1][j] = arrptr[i][B*j];
			}
		}
	}
	return arrptr;
}

template <typename T>
void Tfree_genarray(T  ***arrptr)
{
	if(arrptr == NULL)
	{
		printf("Wrong Input pointer : Tfree_genarray!\n");
		exit(-1);
	}
	T ** ptr_tmp = *arrptr;
	T current_num = ptr_tmp[0][0];
	for(T i = 0; i < current_num; i++)
	{
		delete [] ptr_tmp[i];
	}
	delete [] ptr_tmp;
	ptr_tmp = NULL;
	*arrptr = ptr_tmp;
}

template <typename T, typename TF>
T Tfind_arrindexN(T  **arrptr, TF* const kmer_arr,TF *kmer_inf,uint32_t len)
{
	if(arrptr == NULL)
	{
		printf("Wrong Input : Tfind_arrindexN!\n");
	}
	T level_num = arrptr[0][0]-1;
	T current_num;
	TF *pkmer_cur;
	T lowbound = 1,upbound = B;
	for(T i = 1; i <= level_num; i++)
	{
		current_num = arrptr[i][0];
		if(upbound > current_num)
		{
			upbound = current_num;
		}
		for(T j = lowbound; j <= upbound; j++)
		{
			pkmer_cur = kmer_arr + len*(arrptr[i][j]-1);
			if(cmp256BitKmer(pkmer_cur,kmer_inf,len) != 0)
			{
				lowbound = (j - 1)*B + 1;
				upbound = B*j;
				break;
			}
		}
	}
	lowbound--,upbound--;
	for(T k = lowbound; k <= upbound; k++)
	{
		TF * pkmer_tmp = kmer_arr+len*k;
		if(cmp256BitKmer(pkmer_tmp,kmer_inf,len) == 2)
		{
			return k;
		}
		if(cmp256BitKmer(pkmer_tmp,kmer_inf,len) == 1)
		{
			break;
		}
	}
	return ULLONG_MAX;

}

#endif /* BINARY_SEARCH_H_ */
