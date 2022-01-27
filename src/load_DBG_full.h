/*
 * load_DBG_full.h
 *
 *  Created on: Oct 9, 2019
 *      Author: bio
 */

#ifndef LOAD_DBG_FULL_H_
#define LOAD_DBG_FULL_H_

#include "basic.h"

struct unipath
{
	uint32_t start;
	uint32_t len;
	uint8_t ref_id;
};

struct kmer_detail
{
	uint8_t is_branched;
	uint8_t ad;
	uint32_t unipath_id;
	uint16_t unipath_offset;
};

struct dBG
{
	struct NodeBit** p_kmer_root;
	struct kmer_detail *p_kmer_detail;
	char ** p_unipath;
	uint32_t kN;
	uint32_t uN;
	uint32_t L;

};

struct para_merge
{
	uint64_t *p_start_kmer;
	uint64_t *p_des_kmer;
	uint32_t *p_start_id;
	uint32_t *p_des_id;
	uint8_t *p_start_ad;
	uint8_t *p_des_ad;
	uint32_t blocksize;
	uint32_t thread_num;
	uint32_t thread_id;
	uint64_t ukN;
	uint32_t unitperkmer;
	uint32_t N_block;

};

struct para_dBGindex
{
	uint64_t bkN;
	uint64_t *p_branchedkmer;
	uint8_t *p_branchedkmerad;
	char **upath_arr;
	uint32_t uN;
	uint32_t upathid_start0;
	uint64_t ukN;
	uint64_t *p_unbranchedkmer;
	uint32_t *p_unbranchedkmerid;
	uint64_t **p2_bkmer;
	uint64_t **p2_ukmer;
};

extern struct para_dBGindex gdBGindex;

void get_unipath_struct(uint64_t a,struct unipath *b);
uint32_t get_Dbg_file_name(const char* p_dbg_path, char *** p_dbg_file);
uint64_t get_total_kmer_number(char ** p_dbg_file, uint64_t *kN, uint64_t *ukN);
uint32_t get_total_unipath_number(char ** p_dbg_file);
void get_unipath_kmer_ad(char* p_cur,uint32_t kmer_len,uint8_t * ad);
void save_kmer_details(struct kmer_detail * a,\
		uint8_t is_branched,\
		uint8_t ad,\
		uint32_t unipath_id,\
		uint16_t unipath_offset);

void SortA_umers(struct dBG * p_dBG,char * p_dbg_path,uint32_t thread_num);
void SortFile_umers(struct dBG * p_dBG,char * p_file_path,uint32_t thread_num);
void Gen_navigatSeq(struct dBG * p_dBG, char * p_dbg_path);
//void gen_dBG_index(struct bit256KmerPara bit_para, struct para_dBGindex * p_sdBGidx, char * p_dbg_path, uint32_t thread_num);
void gen_dBG_index(struct bit256KmerPara &bit_para, struct para_dBGindex &sdBGidx, const char * p_dbg_path, uint32_t thread_num, char **reftmp);
void free_dBGindex(struct para_dBGindex &sdBGidx);
void SortFile_umers_forB(char * p_file_path,uint32_t thread_num,uint8_t unitperkmer);
uint64_t * SortFile_umers_for_split(uint64_t unitnum,uint64_t *kmer_array,uint32_t thread_num,uint8_t unitperkmer);

#endif /* LOAD_DBG_FULL_H_ */
