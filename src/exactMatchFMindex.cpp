/*
 * exactMatchFMindex.cpp
 *
 *  Created on: Apr 15, 2021
 *      Author: bio
 */

#include "exactMatchFMindex.h"

uint32_t calc_OCC(const sFMindex &mem,char c,size_t pos)
{
	//calc the occ_array  input : pos(pos in bseq) gap(every gap pos save a occa)
	uint8_t offset = char_to_uint8_table[c];
	if(pos % mem.occ_gap == 0)
	{
		return *(mem.occa + (pos/mem.occ_gap)*4 + offset);
	}
	else
	{
		uint32_t rslt = *(mem.occa + (pos/mem.occ_gap)*4 + offset);
		uint32_t len = pos > mem.occ_gap ? pos - pos % mem.occ_gap : 0;
		for(uint32_t i = 0; i < pos % mem.occ_gap; i++)
		{
			if(mem.b[len+i] == c)
			{
				rslt++;
			}
		}
		return rslt;
	}
}

uint32_t* calc_SArangeSeq(const sFMindex &mem,const char *read)
{
	char ch;
	int32_t i = strlen(read) - 1;
	uint32_t *p = (uint32_t *)malloc(sizeof(uint32_t)*2);
	if(0 > i)
	{
		p[0] = -1;
		p[1] = 0;
		return p;
	}
	ch = read[i];
	p[0] = calc_C(mem,ch) + 1;
	p[1] = calc_C(mem,C_next(ch)) + 1;
	i--;
	while(p[0] <= p[1] && i >= 0)
	{
		ch = read[i];
		p[0] = char_to_uint8_table[ch] != 0x04 ? LF_Mapping_l(mem,ch,p[0]) : -1;
		p[1] = LF_Mapping_h(mem,ch,p[1]);
		i--;
	}
	return p;
}

uint32_t* calc_SArangeChar(const sFMindex &mem,uint32_t *pre, char ch)
{
	uint32_t *p = (uint32_t *)malloc(sizeof(uint32_t)*2);
	p[0] = LF_Mapping_l(mem,ch,pre[0]);
	p[1] = LF_Mapping_h(mem,ch,pre[1]);
	return p;
}

uint32_t calc_SA(const sFMindex &mem,size_t pos)
{
	if(pos % mem.sa_gap == 0)
	{
		return mem.sa[pos/mem.sa_gap];
	}
	else
	{
		size_t tmp = pos;
		uint32_t cnt = 0;
		char ch = mem.b[tmp];
		while(tmp % mem.sa_gap != 0)
		{
			if(ch == '$')
			{
				return cnt;
			}
			tmp = LF_Mapping(mem,ch,tmp);
			ch = mem.b[tmp];
			cnt++;
		}
		return mem.sa[tmp/mem.sa_gap] + cnt;
	}
}




