/*
 * BPlusTree_full.cpp
 *
 *  Created on: Sep 25, 2019
 *      Author: pluto
 */
#include "Binary_Search.h"

uint32_t  **generate_array(uint32_t total_num, uint32_t save2file)  // save2file == 1  save arr to file
{
	if(total_num <= 0 || total_num < B)
	{
		printf("Wrong Input Number : generate_array!\n");
		exit(-1);
	}
	uint32_t level_num;
	float rslt;
	rslt = log10(total_num) / log10(B);
	level_num = (uint32_t)ceil(rslt);
	uint32_t ** arrptr;
	arrptr = (uint32_t **)malloc(level_num * sizeof(uint32_t *));
	arrptr[0] = (uint32_t *)malloc(1 * sizeof(uint32_t));
	arrptr[0][0] = level_num;
	uint32_t remain_num,current_num;
	current_num = total_num;
	for(uint32_t i = level_num; i > 1; i--)
	{
		if(i == level_num)
		{
			remain_num = current_num % B;
		}
		current_num = (uint32_t)ceil((double)current_num / B);
		arrptr[i-1] = (uint32_t *)malloc((current_num+1) * sizeof(uint32_t));
		arrptr[i-1][0] = current_num;
		for(uint32_t j = 1; j <= current_num; j++)
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
	if(save2file)
	{
		FILE * fp = NULL;
		uint32_t save_len = 0;
		fp = fopen("./BranchArray","wb+");
		fwrite(&arrptr[0][0],sizeof(uint32_t),1,fp);
		for(uint32_t i = 1; i < level_num; i++)
		{
			save_len = arrptr[i][0];
			fwrite(&arrptr[i][0],sizeof(uint32_t),save_len+1,fp);
		}
		fclose(fp);
	}
	return arrptr;
}

void free_genarray(uint32_t  ***arrptr)
{
	if(arrptr == NULL)
	{
		printf("Wrong Input pointer : free_array!\n");
		exit(-1);
	}
	uint32_t ** ptr_tmp = *arrptr;
	uint32_t current_num = ptr_tmp[0][0];
	for(uint32_t i = 0; i < current_num; i++)
	{
		free(ptr_tmp[i]);
	}
	free(ptr_tmp);
	ptr_tmp = NULL;
	*arrptr = ptr_tmp;
}

uint32_t find_arrindex(uint32_t  **arrptr, uint64_t* const kmer_arr,uint64_t kmer_inf)
{
	if(arrptr == NULL)
	{
		printf("Wrong Input : BPlusTree_full Array!\n");
	}
	uint32_t level_num = arrptr[0][0]-1;
	uint32_t current_num;
	uint64_t kmer_cur;
	uint32_t lowbound = 1,upbound = B;
	for(uint32_t i = 1; i <= level_num; i++)
	{
		current_num = arrptr[i][0];
		if(upbound > current_num)
		{
			upbound = current_num;
		}
		for(uint32_t j = lowbound; j <= upbound; j++)
		{
			kmer_cur = kmer_arr[arrptr[i][j]-1];
			if(kmer_cur >= kmer_inf )
			{
				lowbound = (j - 1)*B + 1;
				upbound = B*j;
				break;
			}
		}
	}
	lowbound--,upbound--;
	for(uint32_t k = lowbound; k <= upbound; k++)
	{
		if(kmer_arr[k] == kmer_inf)
		{
			return k;
		}
		if(kmer_arr[k] > kmer_inf)
		{
			break;
		}
	}
	return UINT_MAX;

}

uint32_t find_arrindexN(uint32_t  **arrptr, uint64_t* const kmer_arr,uint64_t *kmer_inf,uint32_t len)
{
	if(arrptr == NULL)
	{
		printf("Wrong Input : BPlusTree_full Array!\n");
	}
	uint32_t level_num = arrptr[0][0]-1;
	uint32_t current_num;
	uint64_t *pkmer_cur;
	uint32_t lowbound = 1,upbound = B;
	for(uint32_t i = 1; i <= level_num; i++)
	{
		current_num = arrptr[i][0];
		if(upbound > current_num)
		{
			upbound = current_num;
		}
		for(uint32_t j = lowbound; j <= upbound; j++)
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
	for(uint32_t k = lowbound; k <= upbound; k++)
	{
		uint64_t * pkmer_tmp = kmer_arr+len*k;
		if(cmp256BitKmer(pkmer_tmp,kmer_inf,len) == 2)
		{
			return k;
		}
		if(cmp256BitKmer(pkmer_tmp,kmer_inf,len) == 1)
		{
			break;
		}
	}
	return UINT_MAX;

}

uint32_t find_supathindex(uint32_t  **arrptr, uint32_t *super_arr,uint32_t pos)
{
	if(arrptr == NULL || super_arr == NULL)
	{
		printf("Wrong Input : find_supathindex!\n");
	}
	uint32_t level_num = arrptr[0][0]-1;
	uint32_t current_num;
	uint64_t super_cur;
	uint32_t lowbound = 1,upbound = B;
	uint32_t index_tmp = 0;
	for(uint32_t i = 1; i <= level_num; i++)
	{
		current_num = arrptr[i][0];
		if(upbound > current_num)
		{
			upbound = current_num;
		}
		for(uint32_t j = lowbound; j <= upbound; j++)
		{
			super_cur = super_arr[arrptr[i][j]-1];
			if(super_cur >= pos)
			{
				lowbound = (j - 1)*B + 1;
				upbound = B*j;
				break;
			}
		}
	}
	lowbound--,upbound--;
	for(uint32_t k = lowbound; k <= upbound; k++)
	{
		if(super_arr[k] < pos)
		{
			index_tmp = k+1;
		}
	}
	return index_tmp + 1;
}

