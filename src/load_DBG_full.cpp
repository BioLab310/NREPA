/*
 * load_DBG_full.cpp
 *
 *  Created on: Oct 9, 2019
 *      Author: bio
 */

#include "load_DBG_full.h"
#include "Binary_Search.h"
#include "BplusTreeBit.h"

struct para_dBGindex gdBGindex;

void get_unipath_struct(uint64_t a,struct unipath *b)
{
	b->start = a >> 25;
	b->len = (a >> 8) & 0x1ffff;
	b->ref_id = a & 0xff;
}

uint32_t get_Dbg_file_name(const char* p_dbg_path, char *** p_dbg_file)
{
	uint32_t fN = 0;
	FILE * fp = NULL; //文件指针。
    int32_t c = 0, lc = 0; //c为文件当前字符，lc为上一个字符，供结尾判断用。
    uint32_t line_num = 0; //行数统计
    fp = fopen(p_dbg_path, "r");//以只读方式打开文件。
    if(fp == NULL)
    {
    	printf("%s is not exist!\n",p_dbg_path);
    	return 0;
    }
    while((c = fgetc(fp)) != EOF) //逐个读入字符直到文件结尾
    {
        if(c == '\n') //统计行数。
        {
        	line_num++;
        }
        lc = c; //保存上一字符。
    }
    if(lc != '\n')//处理末行
    {
    	line_num++;
    }
    rewind(fp); //重新指向文件开始
    char ** p_tmp = (char **)malloc(sizeof(char *) * line_num);
//    char path_tmp[32] = {0};
    for( ; fN < line_num; fN++)
    {
    	char *name_tmp = (char *)malloc(sizeof(char) * 32);
		fgets(name_tmp,32,fp);
//		strcat(name_tmp,path_tmp);
		name_tmp[strlen(name_tmp)-1] = '\0'; // replace '\n' with '\0'
		*(p_tmp+fN) = name_tmp;
    }
    fclose(fp);
    *p_dbg_file = p_tmp;
	return fN;
}

uint64_t get_total_kmer_number(char ** p_dbg_file, uint64_t *kN, uint64_t *ukN)
{
	uint64_t tN = 0;
	uint64_t kN_b = 0;
	FILE * fp0 = NULL;
	FILE * fp2 = NULL;
    uint64_t len = 0; //c为文件当前字符，lc为上一个字符，供结尾判断用。
	fp0 = fopen(p_dbg_file[0], "rb");
	fp2 = fopen(p_dbg_file[2], "rb");
	if(fp0 == NULL || fp2 == NULL)
	{
		printf("%s or %s not exist!\n",p_dbg_file[0],p_dbg_file[2]);
		return 0;
	}
	fseek(fp0,0,SEEK_END);//将文件内部的指针指向文件末尾
	len = ftell(fp0);//获取文件长度，（得到文件位置指针当前位置相对于文件首的偏移字节数）
	kN_b = len / sizeof(uint64_t);
    fclose(fp0);

    fseek(fp2,0,SEEK_END);//将文件内部的指针指向文件末尾
    len = ftell(fp2);
    fseek(fp2,0,0);
//    cout << len << endl;
    uint64_t *ptmp = NULL;
    ptmp = (uint64_t *)malloc(len);
    memset(ptmp,0,len);
    fread(ptmp,1,len,fp2);
    len = len / sizeof(uint64_t);
    uint32_t kN_unipath;
    kN_unipath=0;
    for(uint32_t i = 0; i < len; i++)
    {
    	uint64_t tmp;
    	tmp=(*(ptmp+i) >> 8 ) & 0x1ffff;
    	kN_unipath = kN_unipath + tmp;
    }

    tN = kN_b+kN_unipath;
    *kN = kN_b;
    *ukN = kN_unipath;

    free(ptmp);
    fclose(fp2);
	return tN;
}

uint32_t get_total_unipath_number(char ** p_dbg_file)
{
	uint32_t uN = 0;
	FILE * fp = NULL;
    uint32_t len = 0;
	fp = fopen(p_dbg_file[2],"r");
	if(fp == NULL)
	{
		printf("%s is not exist!\n",p_dbg_file[2]);
		return 0;
	}
	fseek(fp,0,SEEK_END);//将文件内部的指针指向文件末尾
	len = ftell(fp);//获取文件长度，（得到文件位置指针当前位置相对于文件首的偏移字节数）
	uN = len / sizeof(uint64_t);
    fclose(fp);
	return uN;
}

void get_unipath_kmer_ad(char* p_cur,uint32_t kmer_len,uint8_t * ad)
{
	*ad = 0;
	switch(*(p_cur-1))
	{
		case 'A':
			*ad += 0x80;
			break;
		case 'C':
			*ad += 0x40;
			break;
		case 'G':
			*ad += 0x20;
			break;
		case 'T':
			*ad += 0x10;
			break;
		default:
			break;
	}
	switch(*(p_cur+kmer_len))
	{
		case 'A':
			*ad += 0x08;
			break;
		case 'C':
			*ad += 0x04;
			break;
		case 'G':
			*ad += 0x02;
			break;
		case 'T':
			*ad += 0x01;
			break;
		default:
			break;
	}
}

void save_kmer_details(struct kmer_detail * a,\
		uint8_t is_branched,\
		uint8_t ad,\
		uint32_t unipath_id,\
		uint16_t unipath_offset)
{
	a->is_branched = is_branched;
	a->ad = ad;
	a->unipath_id = unipath_id;
	a->unipath_offset = unipath_offset;
}

void *merge(void *arg)
{
	struct para_merge* tmp=(struct para_merge*)arg;

	uint32_t groupsize;
	groupsize=2*tmp->thread_num;
	uint32_t N_group;
	N_group=ceil((double)tmp->N_block/groupsize);
	uint32_t threadsize;
	threadsize=2;
	for(uint32_t i=0;i<N_group;i++)
	{
		uint64_t start;
		start=i*groupsize+tmp->thread_id*threadsize;

		if(start>=tmp->N_block)
		{
			return NULL;
		}
		else if(start==tmp->N_block-1)
		{
			//copy
			uint32_t offset=tmp->blocksize*(tmp->N_block-1);
			uint32_t last_kmer=tmp->ukN-offset;
			uint64_t *p_des_kmer_tmp = tmp->p_des_kmer + offset*tmp->unitperkmer;
			uint64_t *p_start_kmer_tmp = tmp->p_start_kmer + offset*tmp->unitperkmer;
			uint32_t *p_des_id_tmp = tmp->p_des_id + offset;
			uint32_t *p_start_id_tmp = tmp->p_start_id + offset;
			uint8_t *p_des_ad_tmp = tmp->p_des_ad + offset;
			uint8_t *p_start_ad_tmp = tmp->p_start_ad + offset;
			memcpy(p_des_kmer_tmp,p_start_kmer_tmp,last_kmer*tmp->unitperkmer*sizeof(uint64_t));
			if(tmp->p_start_id != NULL && tmp->p_des_id != NULL)
			{
				memcpy(p_des_id_tmp,p_start_id_tmp,last_kmer*sizeof(uint32_t));
			}
			if(tmp->p_start_ad != NULL && tmp->p_des_ad != NULL)
			{
				memcpy(p_des_ad_tmp,p_start_ad_tmp,last_kmer*sizeof(uint8_t));
			}

		}
		else
		{
			//merge
			uint64_t *p1,*p2,*p3;
			uint32_t *id1,*id2,*id3;

			uint8_t *ad1,*ad2,*ad3;

			p1=tmp->p_start_kmer+start*tmp->blocksize*tmp->unitperkmer;
			p2=tmp->p_start_kmer+(start+1)*tmp->blocksize*tmp->unitperkmer;
			p3=tmp->p_des_kmer+start*tmp->blocksize*tmp->unitperkmer;

			id1=tmp->p_start_id+start*tmp->blocksize;
			id2=tmp->p_start_id+(start+1)*tmp->blocksize;
			id3=tmp->p_des_id+start*tmp->blocksize;

			ad1=tmp->p_start_ad+start*tmp->blocksize;
			ad2=tmp->p_start_ad+(start+1)*tmp->blocksize;
			ad3=tmp->p_des_ad+start*tmp->blocksize;

			uint32_t Size_block_loop_p2;

			if(i == (N_group - 1))
			{
				uint64_t last_p2=(start+1)*tmp->blocksize;
				if(last_p2+tmp->blocksize>tmp->ukN)
				{
					Size_block_loop_p2=tmp->ukN-last_p2;
				}
				else
				{
					Size_block_loop_p2=tmp->blocksize;
				}
			}
			else
			{
				Size_block_loop_p2=tmp->blocksize;
			}

			uint64_t cnt1 = 0, cnt2 = 0, memcp_n = 0;
			while(cnt1 < tmp->blocksize && cnt2 < Size_block_loop_p2)
			{
				if(cmp256BitKmer(p1,p2,tmp->unitperkmer)==0)
				{
					cnt1++;
					kmercpy(p3,p1,tmp->unitperkmer);
					p1+=tmp->unitperkmer;
					p3+=tmp->unitperkmer;

					if(tmp->p_start_id != NULL && tmp->p_des_id != NULL)
					{
						*id3=*id1;
						id1++;
						id3++;
					}

					if(tmp->p_start_ad != NULL && tmp->p_des_ad != NULL)
					{
						*ad3=*ad1;
						ad1++;
						ad3++;
					}
				}
				else
				{
					cnt2++;
					kmercpy(p3,p2,tmp->unitperkmer);
					p2+=tmp->unitperkmer;
					p3+=tmp->unitperkmer;

					if(tmp->p_start_id != NULL && tmp->p_des_id != NULL)
					{
						*id3=*id2;
						id2++;
						id3++;
					}

					if(tmp->p_start_ad != NULL && tmp->p_des_ad != NULL)
					{
						*ad3=*ad2;
						ad2++;
						ad3++;
					}
				}
			}
			if(cnt2==Size_block_loop_p2)
			{
				memcp_n = tmp->blocksize - cnt1;
				memcpy(p3,p1,memcp_n*tmp->unitperkmer*sizeof(uint64_t));
				if(tmp->p_start_id != NULL && tmp->p_des_id != NULL)
				{
					memcpy(id3,id1,memcp_n*sizeof(uint32_t));
				}

				if(tmp->p_start_ad != NULL && tmp->p_des_ad != NULL)
				{
					memcpy(ad3,ad1,memcp_n*sizeof(uint8_t));
				}
			}
			else
			{
				memcp_n = Size_block_loop_p2 - cnt2 ;
				memcpy(p3,p2,memcp_n*tmp->unitperkmer*sizeof(uint64_t));
				if(tmp->p_start_id != NULL && tmp->p_des_id != NULL)
				{
					memcpy(id3,id2,memcp_n*sizeof(uint32_t));
				}

				if(tmp->p_start_ad != NULL && tmp->p_des_ad != NULL)
				{
					memcpy(ad3,ad2,memcp_n*sizeof(uint8_t));
				}
			}
		}
	}


	return NULL;
}

void SortA_umers(struct dBG * p_dBG,char * p_dbg_path,uint32_t thread_num)
{
	struct bit256KmerPara bit_para;
	get_para(&bit_para,p_dBG->L);
	uint32_t unitperkmer;
	unitperkmer = (uint32_t)ceil((double)p_dBG->L / 32);
	cout << "unitperkmer : " << unitperkmer << endl;

	//1)解析文件名：从p_dbg_path文件中读入4个文件名，分别为：
	//kmer026, kmer026ad, unipath026, nomLeu2.fa
	char ** p_dbg_file;

	uint32_t fN;
	fN = get_Dbg_file_name(p_dbg_path,&p_dbg_file);

	uint32_t uN;
	uN = get_total_unipath_number(p_dbg_file);
	char * p_ref;
	uint32_t ref_len;
	uint32_t cur_ref_id;
	struct unipath cur_struct_unipath_tmp;

	//2)读unipath
	FILE * fp_unipath;
	fp_unipath = fopen(p_dbg_file[2],"rb");
	cur_ref_id = 0;
	ReadSeq(&p_ref,&ref_len,p_dbg_file[cur_ref_id+3]);
	uint64_t *p_unipath;
	uint64_t ukN = 0;
	p_unipath=(uint64_t*)malloc(sizeof(uint64_t)*uN);
	fread(p_unipath,sizeof(uint64_t),uN,fp_unipath);
	fclose(fp_unipath);

	//3)计算umer数量并开辟空间
	for(uint32_t i=0;i<uN;i++)
	{
		get_unipath_struct(p_unipath[i],&cur_struct_unipath_tmp);
		if(cur_struct_unipath_tmp.start == 0)
		{
			continue;
		}
		ukN += cur_struct_unipath_tmp.len;
	}
	uint64_t * const hashvalue_tmp = (uint64_t *)malloc(sizeof(uint64_t) * unitperkmer);
	uint64_t * ukmer_array;
	uint64_t unitnum = ukN * unitperkmer;
	ukmer_array = (uint64_t *)malloc(sizeof(uint64_t)*unitnum);
	memset(ukmer_array,0,sizeof(uint64_t)*unitnum);
	uint32_t * ukmerid_array;
	ukmerid_array = (uint32_t *)malloc(sizeof(uint32_t)*ukN);
	memset(ukmerid_array,0,sizeof(uint32_t)*ukN);


	struct timeval tvs,tve;
	gettimeofday(&tvs,NULL);
	//4)分块排序
	//4.2)分块排序
//	char outputunipathName[32] = {0};
//	sprintf(outputunipathName,"%s%s",p_dbg_file[2],"_Seq");
//	FILE * fp_unipathseq;
//	fp_unipathseq = fopen(outputunipathName,"w");
	char * p_unipath_tmp = NULL;
	uint64_t uai = 0, uiai = 0;
	uint32_t start0_unipath_id = -1;
	uint32_t blocksize = (uint32_t)pow(2,1);
	for(uint32_t i=0;i<uN;i++)
	{
		get_unipath_struct(p_unipath[i],&cur_struct_unipath_tmp);
		if(cur_struct_unipath_tmp.ref_id!=cur_ref_id)
		{
			free(p_ref);
			cur_ref_id = cur_struct_unipath_tmp.ref_id;
			ReadSeq(&p_ref,&ref_len,p_dbg_file[cur_ref_id+3]);
		}
		uint32_t cur_unipath_knum = cur_struct_unipath_tmp.len;
		uint32_t cur_unipath_len = (cur_struct_unipath_tmp.len + p_dBG->L-1);
		char * p_char_tmp;
		p_char_tmp = p_ref+cur_struct_unipath_tmp.start;

//		p_unipath_tmp = (char *)malloc(sizeof(char)*cur_unipath_len+3);
//		memset(p_unipath_tmp,'\0',cur_unipath_len);
//		strncpy(p_unipath_tmp,p_char_tmp-1,cur_unipath_len+2);
//		fprintf(fp_unipathseq, "%s\n", p_unipath_tmp);
//		free(p_unipath_tmp);
//		p_unipath_tmp = NULL;
		if(cur_struct_unipath_tmp.start == 0)
		{
			start0_unipath_id = i;
			continue;
		}
		cal_hash_value_directly_256bit(p_char_tmp,hashvalue_tmp,bit_para);
		kmercpy(ukmer_array+uai,hashvalue_tmp,unitperkmer);
		uai=uai+unitperkmer;
		ukmerid_array[uiai++] = i;

		for(uint32_t j=1; j < cur_unipath_knum; j++)
		{
			cal_hash_value_indirectly_256bit(p_char_tmp+j,hashvalue_tmp,hashvalue_tmp,bit_para);
			kmercpy(ukmer_array+uai,hashvalue_tmp,unitperkmer);
			uai=uai+unitperkmer;
			ukmerid_array[uiai++] = i;
		}
	}
	cout << "start0_unipath_id:" << start0_unipath_id << endl;
//	fclose(fp_unipathseq);

	//两两排序
	uint64_t * cmp_p1,*cmp_p2;
	uint32_t *cmp_id1,*cmp_id2;
	uint32_t cmp_idtmp;
	cmp_p1 = ukmer_array;
	cmp_p2 = ukmer_array + unitperkmer;
	cmp_id1 = ukmerid_array;
	cmp_id2 = cmp_id1 + 1;
	while(cmp_p1 != ukmer_array + ukN * unitperkmer && cmp_p2 != ukmer_array + ukN * unitperkmer)
	{
		if(cmp256BitKmer(cmp_p1,cmp_p2,unitperkmer)==1)
		{
			kmercpy(hashvalue_tmp,cmp_p1,unitperkmer);
			kmercpy(cmp_p1,cmp_p2,unitperkmer);
			kmercpy(cmp_p2,hashvalue_tmp,unitperkmer);
			cmp_idtmp = *cmp_id1;
			*cmp_id1 = *cmp_id2;
			*cmp_id2 = cmp_idtmp;
		}
		cmp_p1 += 2*unitperkmer;
		cmp_p2 += 2*unitperkmer;
		cmp_id1 += 2;
		cmp_id2 += 2;
	}



	cout << "uai:" << uai << " uiai:" << uiai << " ukN:" << ukN << endl;
	gettimeofday(&tve,NULL);
	double span = tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;
	cout << "---2exchange sort time is: "<< span << endl;


	gettimeofday(&tvs,NULL);
	//5)整体排序
	//5.1分配轮换空间
	uint64_t *ukmer_array_1;
	ukmer_array_1=(uint64_t *)malloc(sizeof(uint64_t)*unitnum);
	memset(ukmer_array_1,0,sizeof(uint64_t)*unitnum);
	uint32_t * ukmerid_array_1;
	ukmerid_array_1 = (uint32_t *)malloc(sizeof(uint32_t)*ukN);
	memset(ukmerid_array_1,0,sizeof(uint32_t)*ukN);

	//5.2计算合并次数
	uint32_t N_merge;
	uint32_t N_block;
	N_block=ceil((double)ukN/blocksize);
	N_merge=ceil(log(N_block)/log(2));

	uint32_t N_block_tmp=N_block;
	uint64_t *p_start_kmer,*p_des_kmer,*p_kmer_tmp;
	uint32_t *p_start_id,*p_des_id,*p_id_tmp;
	p_start_kmer=ukmer_array;
	p_des_kmer=ukmer_array_1;
	p_start_id=ukmerid_array;
	p_des_id=ukmerid_array_1;
	uint32_t last_kmer;
	for(uint32_t i=0;i<N_merge;i++)
	{
		uint32_t Size_block_loop=blocksize*pow(2,i);
		uint32_t Size_block_loop_p2;

		if(thread_num==1)
		{
			for(uint32_t j=0;j<N_block_tmp/2;j++)
			{

				uint64_t *p1,*p2,*p3;
				uint32_t *id1,*id2,*id3;

				p1=p_start_kmer+(2*j)*Size_block_loop*unitperkmer;
				p2=p_start_kmer+(2*j+1)*Size_block_loop*unitperkmer;
				p3=p_des_kmer+(2*j)*Size_block_loop*unitperkmer;

				id1=p_start_id+(2*j)*Size_block_loop;
				id2=p_start_id+(2*j+1)*Size_block_loop;
				id3=p_des_id+(2*j)*Size_block_loop;

				if(j == (N_block_tmp/2 - 1))
				{
					uint64_t last_p2 = (p2-p_start_kmer)/unitperkmer; //(2*j+1)*Size_block_loop
					if(last_p2 + Size_block_loop > ukN)
					{
						Size_block_loop_p2 = ukN-last_p2;
					}
					else
					{
						Size_block_loop_p2=Size_block_loop;
					}
				}
				else
				{
					Size_block_loop_p2=Size_block_loop;
				}

				uint64_t cnt1 = 0, cnt2 = 0, memcp_n = 0;
				while(cnt1 < Size_block_loop && cnt2 < Size_block_loop_p2)
				{
					if(cmp256BitKmer(p1,p2,unitperkmer)==0)
					{
						cnt1++;
						kmercpy(p3,p1,unitperkmer);
						p1+=unitperkmer;
						p3+=unitperkmer;

						*id3=*id1;
						id1++;
						id3++;
					}
					else
					{
						cnt2++;
						kmercpy(p3,p2,unitperkmer);
						p2+=unitperkmer;
						p3+=unitperkmer;

						*id3=*id2;
						id2++;
						id3++;
					}
				}
				if(cnt2==Size_block_loop_p2)
				{
					memcp_n = (Size_block_loop - cnt1) ;
					memcpy(p3,p1,memcp_n*unitperkmer * sizeof(uint64_t));
					memcpy(id3,id1,memcp_n * sizeof(uint32_t));
				}
				else
				{
					memcp_n = (Size_block_loop_p2 - cnt2) ;
					memcpy(p3,p2,memcp_n *unitperkmer *  sizeof(uint64_t));
					memcpy(id3,id2,memcp_n * sizeof(uint32_t));
				}
			}
			if(N_block_tmp%2!=0)
			{
				uint32_t offset = Size_block_loop*(N_block_tmp/2)*2;
				last_kmer = ukN - offset;
				uint64_t *p_des_kmer_tmp = p_des_kmer + offset*unitperkmer;
				uint64_t *p_start_kmer_tmp = p_start_kmer + offset*unitperkmer;
				uint32_t *p_des_id_tmp = p_des_id + offset;
				uint32_t *p_start_id_tmp = p_start_id + offset;
				memcpy(p_des_kmer_tmp,p_start_kmer_tmp,last_kmer*unitperkmer*sizeof(uint64_t));
				memcpy(p_des_id_tmp,p_start_id_tmp,last_kmer*sizeof(uint32_t));
			}
		}
		else
		{
			pthread_t* t;
			t=(pthread_t*)malloc(sizeof(pthread_t)*thread_num);
			struct para_merge * p_para;
			p_para=(struct para_merge *)malloc(sizeof(struct para_merge)*thread_num);
			for(uint32_t i_para=0;i_para<thread_num;i_para++)
			{
				p_para[i_para].p_start_kmer=p_start_kmer;
				p_para[i_para].p_des_kmer=p_des_kmer;
				p_para[i_para].p_start_id=p_start_id;
				p_para[i_para].p_des_id=p_des_id;
				p_para[i_para].blocksize=Size_block_loop;
				p_para[i_para].thread_num=thread_num;
				p_para[i_para].thread_id=i_para;
				p_para[i_para].ukN=ukN;
				p_para[i_para].unitperkmer=unitperkmer;
				p_para[i_para].N_block=N_block_tmp;
				if(pthread_create(t+i_para, NULL, merge, (void*)(p_para+i_para))!=0)
				{
					cout << "error!" << endl;
				}
			}
			for(uint32_t i_para=0;i_para<thread_num;i_para++)
			{
				pthread_join(t[i_para], NULL);
			}
			free(p_para);
			free(t);
		}

		N_block_tmp=ceil((double)N_block_tmp/2);
		p_kmer_tmp=p_des_kmer;
		p_des_kmer=p_start_kmer;
		p_start_kmer=p_kmer_tmp;

		p_id_tmp=p_des_id;
		p_des_id=p_start_id;
		p_start_id=p_id_tmp;
	}

	gettimeofday(&tve,NULL);
	span = tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;
	cout << "---Overall sort time is: "<< span << endl;


	uint64_t *ptr1,*ptr2;
	ptr1 = p_start_kmer;
	ptr2 = p_start_kmer + unitperkmer;
	uint32_t bigger_cnt = 0, equ_cnt = 0;
	while(1)
	{
		if(ptr2 == p_start_kmer + ukN * unitperkmer)
		{
			break;
		}
		if(cmp256BitKmer(ptr1,ptr2,unitperkmer)==0)
		{
			ptr1 += unitperkmer;
			ptr2 += unitperkmer;
		}
		else
		{
			if(cmp256BitKmer(ptr1,ptr2,unitperkmer)==2)
			{
				equ_cnt++;
				cout << "equ : ";
			}
			else
			{
				bigger_cnt++;
				cout << "bigger : ";
				cout << "*ptr2 : " << *ptr2;
			}
			cout << " index : " << ptr1 - p_start_kmer << " ";
			cout << "value : ";
			for(uint32_t k = 0; k < unitperkmer; k++)
			{
				cout << *ptr1;
				ptr1++;
			}
			cout << endl;
			ptr2 += unitperkmer;
//			break;
		}

	}
	char  SortAoutname[32] = {0};
	sprintf(SortAoutname,"%s_%d_t%d","Sortukmer",p_dBG->L,thread_num);
	FILE * fp_unipathseq;
	fp_unipathseq = fopen(SortAoutname,"wb");
	fwrite(ukmer_array,sizeof(uint64_t),unitnum,fp_unipathseq);
	fwrite(ukmerid_array,sizeof(uint32_t),ukN,fp_unipathseq);
	free(ukmer_array_1);
	free(ukmerid_array_1);
	free(ukmer_array);
	free(ukmerid_array);
	free(hashvalue_tmp);
	free(p_unipath);
}

void SortFile_umers(struct dBG * p_dBG,char * p_file_path,uint32_t thread_num)
{
	struct bit256KmerPara bit_para;
	get_para(&bit_para,p_dBG->L);
	uint32_t unitperkmer;
	unitperkmer = (uint32_t)ceil((double)p_dBG->L / 32);
	cout << "unitperkmer : " << unitperkmer << endl;

	uint64_t kN;
	FILE * fp_kmerpath;
	fp_kmerpath = fopen(p_file_path,"rb");
	fseek(fp_kmerpath,0,2);
	kN = ftell(fp_kmerpath)/(unitperkmer*sizeof(uint64_t));
	rewind(fp_kmerpath);

	uint64_t * const hashvalue_tmp = (uint64_t *)malloc(sizeof(uint64_t) * unitperkmer);
	uint64_t * kmer_array;
	uint64_t unitnum = kN * unitperkmer;
	cout << "ByteNum:" << unitnum * sizeof(uint64_t) << endl;
	kmer_array = (uint64_t *)malloc(sizeof(uint64_t)*unitnum);
	memset(kmer_array,0,sizeof(uint64_t)*unitnum);
	fread(kmer_array,sizeof(uint64_t),unitnum,fp_kmerpath);
	//4)分块排序
	//4.2)分块排序
	uint64_t * cmp_p1,*cmp_p2;
	cmp_p1 = kmer_array;
	cmp_p2 = kmer_array + unitperkmer;
	while(cmp_p1 != kmer_array + kN * unitperkmer && cmp_p2 != kmer_array + kN * unitperkmer)
	{
		if(cmp256BitKmer(cmp_p1,cmp_p2,unitperkmer)==1)
		{
			kmercpy(hashvalue_tmp,cmp_p1,unitperkmer);
			kmercpy(cmp_p1,cmp_p2,unitperkmer);
			kmercpy(cmp_p2,hashvalue_tmp,unitperkmer);
		}
		cmp_p1 += 2*unitperkmer;
		cmp_p2 += 2*unitperkmer;
	}

	//5)整体排序
	//5.1分配轮换空间
	uint32_t blocksize = (uint32_t)pow(2,1);
	uint64_t *kmer_array_1;
	kmer_array_1=(uint64_t *)malloc(sizeof(uint64_t)*unitnum);
	memset(kmer_array_1,0,sizeof(uint64_t)*unitnum);

	//5.2计算合并次数
	uint32_t N_merge;
	uint32_t N_block;
	N_block=ceil((double)kN/blocksize);
	N_merge=ceil(log(N_block)/log(2));

	uint32_t N_block_tmp=N_block;
	uint64_t *p_start_kmer,*p_des_kmer,*p_kmer_tmp;
	uint32_t *p_start_id,*p_des_id,*p_id_tmp;
	p_start_kmer = kmer_array;
	p_des_kmer = kmer_array_1;
	p_start_id = NULL;
	p_des_id = NULL;
	uint32_t last_kmer;
	for(uint32_t i=0;i<N_merge;i++)
	{
		uint32_t Size_block_loop=blocksize*pow(2,i);
		uint32_t Size_block_loop_p2;

		if(thread_num==1)
		{
			for(uint32_t j=0;j<N_block_tmp/2;j++)
			{

				uint64_t *p1,*p2,*p3;

				p1=p_start_kmer+(2*j)*Size_block_loop*unitperkmer;
				p2=p_start_kmer+(2*j+1)*Size_block_loop*unitperkmer;
				p3=p_des_kmer+(2*j)*Size_block_loop*unitperkmer;


				if(j == (N_block_tmp/2 - 1))
				{
					uint64_t last_p2 = (p2-p_start_kmer)/unitperkmer;
					if(last_p2 + Size_block_loop > kN)
					{
						Size_block_loop_p2 = kN-last_p2;
					}
					else
					{
						Size_block_loop_p2=Size_block_loop;
					}
				}
				else
				{
					Size_block_loop_p2=Size_block_loop;
				}

				uint64_t cnt1 = 0, cnt2 = 0, memcp_n = 0;
				while(cnt1 < Size_block_loop && cnt2 < Size_block_loop_p2)
				{
					if(cmp256BitKmer(p1,p2,unitperkmer)==0)
					{
						cnt1++;
						kmercpy(p3,p1,unitperkmer);
						p1+=unitperkmer;
						p3+=unitperkmer;
					}
					else
					{
						cnt2++;
						kmercpy(p3,p2,unitperkmer);
						p2+=unitperkmer;
						p3+=unitperkmer;
					}
				}
				if(cnt2==Size_block_loop_p2)
				{
					memcp_n = (Size_block_loop - cnt1) ;
					memcpy(p3,p1,memcp_n*unitperkmer * sizeof(uint64_t));
				}
				else
				{
					memcp_n = (Size_block_loop_p2 - cnt2) ;
					memcpy(p3,p2,memcp_n *unitperkmer *  sizeof(uint64_t));
				}
			}
			if(N_block_tmp%2!=0)
			{
				uint32_t offset = Size_block_loop*(N_block_tmp/2)*2;
				last_kmer = kN - offset;
				uint64_t *p_des_kmer_tmp = p_des_kmer + offset*unitperkmer;
				uint64_t *p_start_kmer_tmp = p_start_kmer + offset*unitperkmer;
				memcpy(p_des_kmer_tmp,p_start_kmer_tmp,last_kmer*unitperkmer*sizeof(uint64_t));
			}
		}
		else
		{
			pthread_t* t;
			t=(pthread_t*)malloc(sizeof(pthread_t)*thread_num);
			struct para_merge * p_para;
			p_para=(struct para_merge *)malloc(sizeof(struct para_merge)*thread_num);
			for(uint32_t i_para=0;i_para<thread_num;i_para++)
			{
				p_para[i_para].p_start_kmer=p_start_kmer;
				p_para[i_para].p_des_kmer=p_des_kmer;
				p_para[i_para].p_start_id=p_start_id;
				p_para[i_para].p_des_id=p_des_id;
				p_para[i_para].blocksize=Size_block_loop;
				p_para[i_para].thread_num=thread_num;
				p_para[i_para].thread_id=i_para;
				p_para[i_para].ukN=kN;
				p_para[i_para].unitperkmer=unitperkmer;
				p_para[i_para].N_block=N_block_tmp;
				if(pthread_create(t+i_para, NULL, merge, (void*)(p_para+i_para))!=0)
				{
					cout << "error!" << endl;
				}
			}
			for(uint32_t i_para=0;i_para<thread_num;i_para++)
			{
				pthread_join(t[i_para], NULL);
			}
			free(p_para);
			free(t);
		}

		N_block_tmp=ceil((double)N_block_tmp/2);
		p_kmer_tmp=p_des_kmer;
		p_des_kmer=p_start_kmer;
		p_start_kmer=p_kmer_tmp;

		p_id_tmp=p_des_id;
		p_des_id=p_start_id;
		p_start_id=p_id_tmp;
	}
	cout << "Sort over !" << endl;

	uint64_t *ptr1,*ptr2;
	ptr1 = p_start_kmer;
	ptr2 = p_start_kmer + unitperkmer;
	uint32_t unequ_cnt = 0;

	while(1)
	{
		if(ptr2 == p_start_kmer + kN * unitperkmer)
		{
			break;
		}
		if(cmp256BitKmer(ptr1,ptr2,unitperkmer)==0)
		{
			ptr1 += unitperkmer;
			ptr2 += unitperkmer;
		}
		else
		{
			unequ_cnt++;
			cout << ptr1 - p_start_kmer << endl;
			ptr1 += unitperkmer;
			ptr2 += unitperkmer;
		}

	}
	cout << "unequ_cnt : " << unequ_cnt << endl;

	char outputFileName[32] = {0};
	sprintf(outputFileName,"%s%s",p_file_path,"_S");
	FILE * fp_output;
	fp_output = fopen(outputFileName,"wb");
	fwrite(p_start_kmer,sizeof(uint64_t)*unitnum,1,fp_output);
	fclose(fp_output);
	free(kmer_array_1);
	free(kmer_array);
	free(hashvalue_tmp);
}

void Gen_navigatSeq(struct dBG * p_dBG, char * p_dbg_path)
{
	//initial
	char ** p_dbg_file;
	get_Dbg_file_name(p_dbg_path,&p_dbg_file);
	struct bit256KmerPara bit_para;
	get_para(&bit_para,p_dBG->L);
	uint32_t unitperkmer;
	unitperkmer = (uint32_t)ceil((double)p_dBG->L / 32);
	cout << "unitperkmer : " << unitperkmer << endl;

	uint64_t kbN;
	FILE* fp_branched_kmer_file;
	fp_branched_kmer_file = fopen(p_dbg_file[0],"rb");
	fseek(fp_branched_kmer_file,0,2);
	kbN = ftell(fp_branched_kmer_file)/(sizeof(uint64_t)*bit_para.kmer64Len);
	fseek(fp_branched_kmer_file,0,0);
	rewind(fp_branched_kmer_file);

	//gen **array char * array  ukmerarray
	uint32_t **bkmer_ptr;
	bkmer_ptr = generate_array(kbN, 0);
	char *p_navigatseq = NULL;
	uint64_t *bkmer_array = NULL;
	bkmer_array = (uint64_t *)malloc(sizeof(uint64_t) * kbN * unitperkmer);
	fread(bkmer_array,sizeof(uint64_t)*unitperkmer,kbN,fp_branched_kmer_file);
	fclose(fp_branched_kmer_file);
	uint64_t * const hashvalue_tmp = (uint64_t *)malloc(sizeof(uint64_t) * unitperkmer);
	uint32_t alloclen = (uint32_t)pow(2,30);
	uint32_t seq_len = alloclen;
	p_navigatseq = (char *)malloc(sizeof(char) * seq_len);
	uint32_t fill_num = 0;
	//read ref seq
	char *p_ref;
	uint32_t ref_len;
	uint32_t cur_ref_id = 0;
	ReadSeq(&p_ref,&ref_len,p_dbg_file[cur_ref_id+3]);
	//traverse ref sequence
	char *p_char_tmp;
	p_char_tmp = p_ref;
	uint32_t offset = 0;
	uint32_t ret;
	cal_hash_value_directly_256bit(p_char_tmp,hashvalue_tmp,bit_para);
	ret = find_arrindex(bkmer_ptr, bkmer_array, *hashvalue_tmp);
	if(ret != UINT_MAX)
	{
		p_navigatseq[fill_num] = *(p_char_tmp+p_dBG->L);
		fill_num++;
	}
	p_char_tmp++;
	offset++;
	while(offset <= ref_len-p_dBG->L)
	{
		cal_hash_value_indirectly_256bit(p_char_tmp,hashvalue_tmp,hashvalue_tmp,bit_para);
		ret = find_arrindex(bkmer_ptr, bkmer_array, *hashvalue_tmp);		//这里有个问题  find_arrindex函数输入uint64_t  L超过32的目前还不可以
		if(ret != UINT_MAX)
		{
			p_navigatseq[fill_num] = *(p_char_tmp+p_dBG->L);
			fill_num++;
		}
		p_char_tmp++;
		offset++;
		if(fill_num == seq_len)
		{
			seq_len += alloclen;
			p_navigatseq = (char *)realloc(p_navigatseq,seq_len);
		}

	}
	printf("fill_num : %d\n",fill_num);
	char outputFileName[32] = {0};
	sprintf(outputFileName,"%s%s",p_dbg_file[0],"_navigatseq");
	FILE *navigatseq_output;
	navigatseq_output = fopen(outputFileName,"wb");
	fwrite(p_navigatseq,sizeof(char),fill_num,navigatseq_output);
	fclose(navigatseq_output);
	free_genarray(&bkmer_ptr);
	free(p_navigatseq);
	free(hashvalue_tmp);
	free(bkmer_array);
}

void Transition(char start,char end,uint8_t* temp_res){
	uint8_t temp_b;
	uint8_t *res;
	res=temp_res;
	switch(start)
	{
		case 'A':
			temp_b=16;
			break;
		case 'C':
			temp_b=32;
			break;
		case 'G':
			temp_b=64;
			break;
		case 'T':
			temp_b=128;
			break;
		case '0':
			temp_b=0;
			break;
		default:
			temp_b=16;
			break;
	}
	(*res)=(*res)|temp_b;
	switch(end)
	{
		case 'A':
			temp_b=1;
			break;
		case 'C':
			temp_b=2;
			break;
		case 'G':
			temp_b=4;
			break;
		case 'T':
			temp_b=8;
			break;
		case '0':
			temp_b=0;
			break;
		default:
			temp_b=1;
			break;
	}
	(*res)=(*res)|temp_b;
}

void gen_dBG_index(struct bit256KmerPara &bit_para, struct para_dBGindex &sdBGidx, const char * p_dbg_path, uint32_t thread_num, char **reftmp)	//用CDBGO生成的文件生成一个整体的dBG index
{
	uint32_t kmerL = bit_para.kmer1Len / 2;
	uint32_t unitperkmer;
	unitperkmer = bit_para.kmer64Len;
//	cout << "unitperkmer : " << unitperkmer << endl;

	//1)解析文件名：从p_dbg_path文件中读入4个文件名，分别为：
		//kmer026, kmer026ad, unipath026, nomLeu2.fa
	char ** p_dbg_file;

	uint32_t fN;
	fN = get_Dbg_file_name(p_dbg_path,&p_dbg_file);

	uint64_t kBN;
	FILE* fp_branched_kmer_file;
	fp_branched_kmer_file = fopen(p_dbg_file[0],"rb");
	fseek(fp_branched_kmer_file,0,2);
	kBN = ftell(fp_branched_kmer_file)/(sizeof(uint64_t)*unitperkmer);
	fseek(fp_branched_kmer_file,0,0);

	uint64_t kBaN;
	FILE* fp_branched_kmerad_file;
	fp_branched_kmerad_file = fopen(p_dbg_file[1],"rb");
	fseek(fp_branched_kmerad_file,0,2);
	kBaN = ftell(fp_branched_kmerad_file)/sizeof(uint8_t);
	fseek(fp_branched_kmerad_file,0,0);
	if(kBN != kBaN)
	{
		cout << p_dbg_file[0] << " kmer number mismatch " << p_dbg_file[1] << endl;
		exit(1);   //return ?
	}
	uint64_t *p_branchedkmer = NULL;
	p_branchedkmer = (uint64_t *)malloc(sizeof(uint64_t) * kBN);
	fread(p_branchedkmer,sizeof(uint64_t)*unitperkmer,kBN,fp_branched_kmer_file);
	fclose(fp_branched_kmer_file);

	uint8_t *p_branchedkmerad = NULL;
	p_branchedkmerad = (uint8_t *)malloc(sizeof(uint8_t) * kBN);
	fread(p_branchedkmerad,sizeof(uint8_t),kBN,fp_branched_kmerad_file);
	fclose(fp_branched_kmerad_file);

	// 1 )save branch kmer inf and branch kmer ad inf to file
	char dBG_index_name[32] = {0};
	char *p_name = p_dbg_file[3];
	char ref_name[16] = {0};
	while(p_name != p_dbg_file[3] + strlen(p_dbg_file[3]))
	{
		if(*p_name == '.')
		{
			break;
		}
		strncat(ref_name,p_name,1);
		p_name++;
	}
	uint64_t ByteCnt = 0;
	sprintf(dBG_index_name,"%s_%d_%s",ref_name,kmerL,"index");
	FILE * fp_dBGindex;
	fp_dBGindex = fopen(dBG_index_name,"wb");
	fwrite(&kBN,sizeof(uint64_t),1,fp_dBGindex);
	sdBGidx.bkN = kBN;
	ByteCnt += sizeof(uint64_t) * 1;
	fwrite(p_branchedkmer,sizeof(uint64_t)*unitperkmer,kBN,fp_dBGindex);
	sdBGidx.p_branchedkmer = p_branchedkmer;
	ByteCnt += sizeof(uint64_t) * unitperkmer * kBN;
	fwrite(p_branchedkmerad,sizeof(uint8_t),kBN,fp_dBGindex);
	sdBGidx.p_branchedkmerad = p_branchedkmerad;
	ByteCnt += sizeof(uint8_t) * kBN;
//	free(p_branchedkmer);	读入内存不能free掉了
	p_branchedkmer = NULL;
//	free(p_branchedkmerad);
	p_branchedkmerad = NULL;
//	cout << "ByteCnt:" << ByteCnt << endl;
//	cout << "branch kmer fwrite done" << endl;

	uint32_t uN;
	uN = get_total_unipath_number(p_dbg_file);
	sdBGidx.uN = uN;
	sdBGidx.upath_arr = (char**)malloc(sizeof(char*)*uN);
	char * p_ref;
	uint32_t ref_len;
	uint32_t cur_ref_id;
	struct unipath cur_struct_unipath_tmp;


	//2)读unipath
	FILE * fp_unipath;
	fp_unipath = fopen(p_dbg_file[2],"rb");
	cur_ref_id = 0;
	ReadSeq(&p_ref,&ref_len,p_dbg_file[cur_ref_id+3]);
	uint64_t *p_unipath;
	uint64_t ukN = 0;
	p_unipath=(uint64_t*)malloc(sizeof(uint64_t)*uN);
	fread(p_unipath,sizeof(uint64_t),uN,fp_unipath);
	fclose(fp_unipath);

	//3)计算umer数量并开辟空间
	for(uint32_t i=0;i<uN;i++)
	{
		get_unipath_struct(p_unipath[i],&cur_struct_unipath_tmp);
//		if(cur_struct_unipath_tmp.start == 0)
//		{
//			continue;
//		}
		ukN += cur_struct_unipath_tmp.len;
	}
	uint64_t * const hashvalue_tmp = (uint64_t *)malloc(sizeof(uint64_t) * unitperkmer);
	uint64_t * ukmer_array;
	uint64_t unitnum = ukN * unitperkmer;
	ukmer_array = (uint64_t *)malloc(sizeof(uint64_t)*unitnum);
	memset(ukmer_array,0,sizeof(uint64_t)*unitnum);
	uint32_t * ukmerid_array;
	ukmerid_array = (uint32_t *)malloc(sizeof(uint32_t)*ukN);
	memset(ukmerid_array,0,sizeof(uint32_t)*ukN);

	//2020.1024 add non-branching kmer adinformation
	uint8_t *ukmerad_array;
	ukmerad_array = (uint8_t *)malloc(sizeof(uint8_t)*ukN);
	memset(ukmerad_array,0,sizeof(uint8_t)*ukN);

	//4)分块排序
	//4.2)分块排序
	char * p_unipath_tmp = NULL;
	uint64_t uai = 0, uiai = 0;
	uint32_t start0_unipath_id = -1;
	uint32_t blocksize = (uint32_t)pow(2,1);
//	cout << "uN:" << uN << " ukN:" << ukN << endl;
	for(uint32_t i=0;i<uN;i++)
	{
		get_unipath_struct(p_unipath[i],&cur_struct_unipath_tmp);
//		if(cur_struct_unipath_tmp.start == 0)
//		{
//			cout << "unipath start from 0" << endl;
//			start0_unipath_id = i;
//			continue;
//		}
		if(cur_struct_unipath_tmp.ref_id!=cur_ref_id)
		{
			free(p_ref);
			cur_ref_id = cur_struct_unipath_tmp.ref_id;
			ReadSeq(&p_ref,&ref_len,p_dbg_file[cur_ref_id+3]);
		}
		uint32_t cur_unipath_knum = cur_struct_unipath_tmp.len;
		uint32_t cur_unipath_len = (cur_struct_unipath_tmp.len + kmerL-1);
		char * p_char_tmp;
		p_char_tmp = p_ref+cur_struct_unipath_tmp.start;
		p_unipath_tmp = (char *)malloc(sizeof(char)*cur_unipath_len+3);
		memset(p_unipath_tmp,'\0',cur_unipath_len+3);
		strncpy(p_unipath_tmp,p_char_tmp-1,cur_unipath_len+2);
		sdBGidx.upath_arr[i] = p_unipath_tmp;
//		fwrite(&cur_unipath_len,sizeof(uint32_t),1,fp_dBGindex);   不需要存  strlen
//		ByteCnt += sizeof(uint32_t) * 1;
		char *cur_upath_ch = p_unipath_tmp;
		fwrite(p_unipath_tmp,sizeof(char),cur_unipath_len+2,fp_dBGindex);
		ByteCnt += sizeof(char) * strlen(p_unipath_tmp);
		cal_hash_value_directly_256bit(p_char_tmp,hashvalue_tmp,bit_para);
		kmercpy(ukmer_array+uai,hashvalue_tmp,unitperkmer);
		uai=uai+unitperkmer;
		Transition(*cur_upath_ch, *(cur_upath_ch + kmerL + 1), ukmerad_array+uiai);	//1024 added
		ukmerid_array[uiai++] = i;
		cur_upath_ch++;		//1024 added

		for(uint32_t j=1; j < cur_unipath_knum; j++)
		{
			cal_hash_value_indirectly_256bit(p_char_tmp+j,hashvalue_tmp,hashvalue_tmp,bit_para);
			kmercpy(ukmer_array+uai,hashvalue_tmp,unitperkmer);
			uai=uai+unitperkmer;
			Transition(*cur_upath_ch, *(cur_upath_ch + kmerL + 1), ukmerad_array+uiai);	//1024 added
			ukmerid_array[uiai++] = i;
			cur_upath_ch++;		//1024 added
		}
		if(cur_upath_ch - p_unipath_tmp != cur_unipath_knum)
		{
			cout << "unipaht i:" << i << endl;
		}
	}
	fwrite(&start0_unipath_id,sizeof(uint32_t),1,fp_dBGindex);
	sdBGidx.upathid_start0 = start0_unipath_id;
	ByteCnt += sizeof(uint32_t) * 1;
	fwrite(&ukN,sizeof(uint64_t),1,fp_dBGindex);
	sdBGidx.ukN = ukN;
	ByteCnt += sizeof(uint64_t) * 1;
	sdBGidx.p2_bkmer = Tgenerate_array(sdBGidx.bkN);
	sdBGidx.p2_ukmer = Tgenerate_array(sdBGidx.ukN);
//	cout << "ByteCnt:" << ByteCnt << endl;
//	cout << "unipath sequence fwrite done" << endl;

		//两两排序
	char ukmer_file[32] = {0};
	char ukmerid_file[32] = {0};
	char ukmerad_file[32] = {0};
	sprintf(ukmer_file,"%s0%d","ukmer",kmerL);
	sprintf(ukmerid_file,"%s0%did","ukmer",kmerL);
	sprintf(ukmerad_file,"%s0%dad","ukmer",kmerL);	//1024 added
	FILE * fp_unbranched_kmer_file = NULL, *fp_unbranched_kmerid_file = NULL, *fp_unbranched_kmerad_file = NULL;
	if(access(ukmer_file,0) != 0 || access(ukmerid_file,0) != 0 || access(ukmerad_file,0) != 0)   //如果ukmer || ukmerid || ukmerad文件不存在  则读入内存并写入文件
	{
		uint64_t * cmp_p1,*cmp_p2;
		uint32_t *cmp_id1,*cmp_id2;
		uint8_t *cmp_ad1,*cmp_ad2;
		uint32_t cmp_idtmp;
		uint8_t cmp_adtmp;
		cmp_p1 = ukmer_array;
		cmp_p2 = ukmer_array + unitperkmer;
		cmp_id1 = ukmerid_array;
		cmp_id2 = cmp_id1 + 1;
		//1024 added
		cmp_ad1 = ukmerad_array;
		cmp_ad2 = ukmerad_array + 1;
		while(cmp_p1 != ukmer_array + ukN * unitperkmer && cmp_p2 != ukmer_array + ukN * unitperkmer)
		{
			if(cmp256BitKmer(cmp_p1,cmp_p2,unitperkmer)==1)
			{
				kmercpy(hashvalue_tmp,cmp_p1,unitperkmer);
				kmercpy(cmp_p1,cmp_p2,unitperkmer);
				kmercpy(cmp_p2,hashvalue_tmp,unitperkmer);
				cmp_idtmp = *cmp_id1;
				*cmp_id1 = *cmp_id2;
				*cmp_id2 = cmp_idtmp;

				cmp_adtmp = *cmp_ad1;
				*cmp_ad1 = *cmp_ad2;
				*cmp_ad2 = cmp_adtmp;
			}
			cmp_p1 += 2*unitperkmer;
			cmp_p2 += 2*unitperkmer;
			cmp_id1 += 2;
			cmp_id2 += 2;

			//1024 added
			cmp_ad1 += 2;
			cmp_ad2 += 2;
		}

		//5)整体排序
		//5.1分配轮换空间
		uint64_t *ukmer_array_1;
		ukmer_array_1=(uint64_t *)malloc(sizeof(uint64_t)*unitnum);
		memset(ukmer_array_1,0,sizeof(uint64_t)*unitnum);
		uint32_t * ukmerid_array_1;
		ukmerid_array_1 = (uint32_t *)malloc(sizeof(uint32_t)*ukN);
		memset(ukmerid_array_1,0,sizeof(uint32_t)*ukN);

		uint8_t *ukmerad_array_1;
		ukmerad_array_1 = (uint8_t *)malloc(sizeof(uint8_t)*ukN);
		memset(ukmerad_array_1,0,sizeof(uint8_t)*ukN);

		//5.2计算合并次数
		uint32_t N_merge;
		uint32_t N_block;
		N_block=ceil((double)ukN/blocksize);
		N_merge=ceil(log(N_block)/log(2));

		uint32_t N_block_tmp=N_block;
		uint64_t *p_start_kmer,*p_des_kmer,*p_kmer_tmp;
		uint32_t *p_start_id,*p_des_id,*p_id_tmp;

		uint8_t *p_start_ad,*p_des_ad,*p_ad_tmp;	//1024 added

		p_start_kmer = ukmer_array;
		p_des_kmer = ukmer_array_1;
		p_start_id = ukmerid_array;
		p_des_id = ukmerid_array_1;

		p_start_ad = ukmerad_array;				//1024 added
		p_des_ad = ukmerad_array_1;
		uint32_t last_kmer;
		for(uint32_t i=0;i<N_merge;i++)
		{
			uint32_t Size_block_loop=blocksize*pow(2,i);
			uint32_t Size_block_loop_p2;

			if(thread_num==1)
			{
				for(uint32_t j=0;j<N_block_tmp/2;j++)
				{

					uint64_t *p1,*p2,*p3;
					uint32_t *id1,*id2,*id3;
					uint8_t *ad1,*ad2,*ad3;

					p1 = p_start_kmer+(2*j)*Size_block_loop*unitperkmer;
					p2 = p_start_kmer+(2*j+1)*Size_block_loop*unitperkmer;
					p3 = p_des_kmer+(2*j)*Size_block_loop*unitperkmer;

					id1 = p_start_id+(2*j)*Size_block_loop;
					id2 = p_start_id+(2*j+1)*Size_block_loop;
					id3 = p_des_id+(2*j)*Size_block_loop;

					//1024 added
					ad1 = p_start_ad+(2*j)*Size_block_loop;
					ad2 = p_start_ad+(2*j+1)*Size_block_loop;
					ad3 = p_des_ad+(2*j)*Size_block_loop;

					if(j == (N_block_tmp/2 - 1))
					{
						uint64_t last_p2 = (p2-p_start_kmer)/unitperkmer; //(2*j+1)*Size_block_loop
						if(last_p2 + Size_block_loop > ukN)
						{
							Size_block_loop_p2 = ukN-last_p2;
						}
						else
						{
							Size_block_loop_p2=Size_block_loop;
						}
					}
					else
					{
						Size_block_loop_p2=Size_block_loop;
					}

					uint64_t cnt1 = 0, cnt2 = 0, memcp_n = 0;
					while(cnt1 < Size_block_loop && cnt2 < Size_block_loop_p2)
					{
						if(cmp256BitKmer(p1,p2,unitperkmer)==0)
						{
							cnt1++;
							kmercpy(p3,p1,unitperkmer);
							p1+=unitperkmer;
							p3+=unitperkmer;

							*id3=*id1;
							id1++;
							id3++;

							//1024 added
							*ad3=*ad1;
							ad1++;
							ad3++;
						}
						else
						{
							cnt2++;
							kmercpy(p3,p2,unitperkmer);
							p2+=unitperkmer;
							p3+=unitperkmer;

							*id3=*id2;
							id2++;
							id3++;

							//1024 added
							*ad3=*ad2;
							ad2++;
							ad3++;
						}
					}
					if(cnt2==Size_block_loop_p2)
					{
						memcp_n = (Size_block_loop - cnt1);
						memcpy(p3,p1,memcp_n*unitperkmer * sizeof(uint64_t));
						memcpy(id3,id1,memcp_n * sizeof(uint32_t));
						memcpy(ad3,ad1,memcp_n * sizeof(uint8_t));	//1024 added
					}
					else
					{
						memcp_n = (Size_block_loop_p2 - cnt2);
						memcpy(p3,p2,memcp_n *unitperkmer *  sizeof(uint64_t));
						memcpy(id3,id2,memcp_n * sizeof(uint32_t));
						memcpy(ad3,ad2,memcp_n * sizeof(uint8_t));	//1024 added
					}
				}
				if(N_block_tmp%2!=0)
				{
					uint32_t offset = Size_block_loop*(N_block_tmp/2)*2;
					last_kmer = ukN - offset;
					uint64_t *p_des_kmer_tmp = p_des_kmer + offset*unitperkmer;
					uint64_t *p_start_kmer_tmp = p_start_kmer + offset*unitperkmer;
					uint32_t *p_des_id_tmp = p_des_id + offset;
					uint32_t *p_start_id_tmp = p_start_id + offset;

					//1024 added
					uint8_t *p_des_ad_tmp = p_des_ad + offset;
					uint8_t *p_start_ad_tmp = p_start_ad + offset;
					memcpy(p_des_kmer_tmp,p_start_kmer_tmp,last_kmer*unitperkmer*sizeof(uint64_t));
					memcpy(p_des_id_tmp,p_start_id_tmp,last_kmer*sizeof(uint32_t));
					memcpy(p_des_ad_tmp,p_start_ad_tmp,last_kmer*sizeof(uint8_t));
				}
			}
			else
			{
				pthread_t* t;
				t=(pthread_t*)malloc(sizeof(pthread_t)*thread_num);
				struct para_merge * p_para;
				p_para=(struct para_merge *)malloc(sizeof(struct para_merge)*thread_num);
				for(uint32_t i_para=0;i_para<thread_num;i_para++)
				{
					p_para[i_para].p_start_kmer = p_start_kmer;
					p_para[i_para].p_des_kmer = p_des_kmer;
					p_para[i_para].p_start_id = p_start_id;
					p_para[i_para].p_des_id = p_des_id;
					p_para[i_para].p_start_ad = p_start_ad;		//1024 added
					p_para[i_para].p_des_ad = p_des_ad;		//1024 added
					p_para[i_para].blocksize = Size_block_loop;
					p_para[i_para].thread_num = thread_num;
					p_para[i_para].thread_id = i_para;
					p_para[i_para].ukN = ukN;
					p_para[i_para].unitperkmer = unitperkmer;
					p_para[i_para].N_block = N_block_tmp;
					if(pthread_create(t+i_para, NULL, merge, (void*)(p_para+i_para))!=0)
					{
						cout << "error!" << endl;
					}
				}
				for(uint32_t i_para=0;i_para<thread_num;i_para++)
				{
					pthread_join(t[i_para], NULL);
				}
				free(p_para);
				free(t);
			}

			N_block_tmp = ceil((double)N_block_tmp/2);
			p_kmer_tmp = p_des_kmer;
			p_des_kmer = p_start_kmer;
			p_start_kmer = p_kmer_tmp;

			p_id_tmp = p_des_id;
			p_des_id = p_start_id;
			p_start_id = p_id_tmp;

			p_ad_tmp = p_des_ad;
			p_des_ad = p_start_ad;
			p_start_ad = p_ad_tmp;
		}
		fwrite(p_start_kmer,sizeof(uint64_t)*unitperkmer,ukN,fp_dBGindex);
		sdBGidx.p_unbranchedkmer = p_start_kmer;
		fwrite(p_start_id,sizeof(uint32_t),ukN,fp_dBGindex);
		sdBGidx.p_unbranchedkmerid = p_start_id;
		if(p_start_kmer == ukmer_array_1)
		{
			free(ukmer_array);
			free(ukmerid_array);
			free(ukmerad_array);
		}
		else
		{
			free(ukmer_array_1);
			free(ukmerid_array_1);
			free(ukmerad_array_1);
		}
		fp_unbranched_kmer_file = fopen(ukmer_file,"wb");
		fp_unbranched_kmerid_file = fopen(ukmerid_file,"wb");
		fp_unbranched_kmerad_file = fopen(ukmerad_file,"wb");
		fwrite(sdBGidx.p_unbranchedkmer, sizeof(uint64_t), unitnum, fp_unbranched_kmer_file);
		fwrite(sdBGidx.p_unbranchedkmerid, sizeof(uint32_t), ukN, fp_unbranched_kmerid_file);
		fwrite(p_start_ad, sizeof(uint8_t), ukN, fp_unbranched_kmerad_file);
		fclose(fp_unbranched_kmerad_file);
//		cout << "ukmer file size:" << sizeof(uint64_t) * unitnum << "\tukmerid file size:" << sizeof(uint32_t) * ukN << endl;
//		cout << "write ukmer file done !" << endl;
	}
	else
	{
		uint64_t ukNf;
		fp_unbranched_kmer_file = fopen(ukmer_file,"rb");
		fseek(fp_unbranched_kmer_file,0,2);
		ukNf = ftell(fp_unbranched_kmer_file)/(sizeof(uint64_t)*unitperkmer);
		fseek(fp_unbranched_kmer_file,0,0);

		uint64_t ukNfi;
		fp_unbranched_kmerid_file = fopen(ukmerid_file,"rb");
		fseek(fp_unbranched_kmerid_file,0,2);
		ukNfi = ftell(fp_unbranched_kmerid_file)/sizeof(uint32_t);
		fseek(fp_unbranched_kmerid_file,0,0);

		if(ukN != ukNf)
		{
			cout <<"ukN not match ukNf" << endl;
		}
		if(ukNf != ukNfi)
		{
			cout << ukmer_file << " kmer number mismatch " << ukmerid_file << endl;
			exit(1);
		}
		fread(ukmer_array,sizeof(uint64_t)*unitperkmer,ukN,fp_unbranched_kmer_file);
		fread(ukmerid_array,sizeof(uint32_t),ukN,fp_unbranched_kmerid_file);
		sdBGidx.p_unbranchedkmer = ukmer_array;
		sdBGidx.p_unbranchedkmerid = ukmerid_array;
//		cout << "read ukmer from file done !" << endl;
	}
	ByteCnt += sizeof(uint64_t) * unitnum;
	ByteCnt += sizeof(uint32_t) * ukN;
	fclose(fp_unbranched_kmer_file);
	fclose(fp_unbranched_kmerid_file);
//	cout << "ByteCnt:" << ByteCnt << endl;
//	cout << "generate dBG index done" << endl;
//	free(p_ref); //mem-leak !
	for(int i = 0; i < fN; ++i)
	{
		free(p_dbg_file[i]);
		p_dbg_file[i] = NULL;

	}
	free(p_dbg_file);
	p_dbg_file = NULL;
	*reftmp = p_ref;
	free(hashvalue_tmp);
	free(p_unipath);
}

void free_dBGindex(struct para_dBGindex &sdBGidx)
{
	if(sdBGidx.p_branchedkmer != NULL)
	{
		free(sdBGidx.p_branchedkmer);
		sdBGidx.p_branchedkmer = NULL;
		free(sdBGidx.p_branchedkmerad);
		sdBGidx.p_branchedkmerad = NULL;
	}
	if(sdBGidx.p_unbranchedkmer != NULL)
	{
		free(sdBGidx.p_unbranchedkmer);
		sdBGidx.p_unbranchedkmer = NULL;
		free(sdBGidx.p_unbranchedkmerid);
		sdBGidx.p_unbranchedkmerid = NULL;
	}
	if(sdBGidx.upath_arr != NULL)
	{
		for(uint32_t i = 0; i < sdBGidx.uN; i++)
		{
			free(sdBGidx.upath_arr[i]);
			sdBGidx.upath_arr[i] = NULL;
		}
		sdBGidx.upath_arr = NULL;
	}
	if(sdBGidx.p2_bkmer)
	{
		Tfree_genarray(&sdBGidx.p2_bkmer);
	}
	if(sdBGidx.p2_ukmer)
	{
		Tfree_genarray(&sdBGidx.p2_ukmer);
	}
}

void load_dBG_index(struct bit256KmerPara &bit_para, struct para_dBGindex &sdBGidx, const char * p_dbg_path)	//load dBG index 后续工作...
{
	//not used
	uint32_t kmerL = bit_para.kmer1Len / 2;
	uint32_t unitperkmer;
	unitperkmer = bit_para.kmer64Len;

	char ** p_dbg_file;

	uint32_t fN;
	fN = get_Dbg_file_name(p_dbg_path,&p_dbg_file);

	//1)judge whether kBN is equal to kBaN
	uint64_t kBN;
	FILE* fp_branched_kmer_file;
	fp_branched_kmer_file = fopen(p_dbg_file[0],"rb");
	fseek(fp_branched_kmer_file,0,2);
	kBN = ftell(fp_branched_kmer_file)/(sizeof(uint64_t)*unitperkmer);
	fseek(fp_branched_kmer_file,0,0);

	uint64_t kBaN;
	FILE* fp_branched_kmerad_file;
	fp_branched_kmerad_file = fopen(p_dbg_file[1],"rb");
	fseek(fp_branched_kmerad_file,0,2);
	kBaN = ftell(fp_branched_kmerad_file)/sizeof(uint8_t);
	fseek(fp_branched_kmerad_file,0,0);
	if(kBN != kBaN)
	{
		cout << p_dbg_file[0] << " kmer number mismatch " << p_dbg_file[1] << endl;
		exit(1);   //return ?
	}

	uint64_t *p_branchedkmer = NULL;
	p_branchedkmer = (uint64_t *)malloc(sizeof(uint64_t) * kBN);
	fread(p_branchedkmer,sizeof(uint64_t)*unitperkmer,kBN,fp_branched_kmer_file);
	fclose(fp_branched_kmer_file);

	uint8_t *p_branchedkmerad = NULL;
	p_branchedkmerad = (uint8_t *)malloc(sizeof(uint8_t) * kBN);
	fread(p_branchedkmerad,sizeof(uint8_t),kBN,fp_branched_kmerad_file);
	fclose(fp_branched_kmerad_file);

	for(int i = 0; i < fN; ++i)
	{
		free(p_dbg_file[i]);
		p_dbg_file[i] = NULL;
	}
	free(p_dbg_file);
	p_dbg_file = NULL;
}

void SortFile_umers_forB(char * p_file_path,uint32_t thread_num,uint8_t unitperkmer)
{
	uint64_t kN;
	FILE * fp_kmerpath;
	fp_kmerpath = fopen(p_file_path,"rb");
	fseek(fp_kmerpath,0,2);
	kN = ftell(fp_kmerpath)/(unitperkmer*sizeof(uint64_t));
	rewind(fp_kmerpath);

	uint64_t * const hashvalue_tmp = (uint64_t *)malloc(sizeof(uint64_t) * unitperkmer);
	uint64_t * kmer_array;
	uint64_t unitnum = kN * unitperkmer;
	cout << "unitnum:" << unitnum << endl;
	cout << "ByteNum:" << unitnum * sizeof(uint64_t) << endl;
	kmer_array = (uint64_t *)malloc(sizeof(uint64_t)*unitnum);
	memset(kmer_array,0,sizeof(uint64_t)*unitnum);
	fread(kmer_array,sizeof(uint64_t),unitnum,fp_kmerpath);
	//4)分块排序
	//4.2)分块排序
	uint64_t * cmp_p1,*cmp_p2;
	cmp_p1 = kmer_array;
	cmp_p2 = kmer_array + unitperkmer;
	while(cmp_p1 != kmer_array + kN * unitperkmer && cmp_p2 != kmer_array + kN * unitperkmer)
	{
		if(cmp256BitKmer(cmp_p1,cmp_p2,unitperkmer)==1)
		{
			kmercpy(hashvalue_tmp,cmp_p1,unitperkmer);
			kmercpy(cmp_p1,cmp_p2,unitperkmer);
			kmercpy(cmp_p2,hashvalue_tmp,unitperkmer);
		}
		cmp_p1 += 2*unitperkmer;
		cmp_p2 += 2*unitperkmer;
	}

	//5)整体排序
	//5.1分配轮换空间
	uint32_t blocksize = (uint32_t)pow(2,1);
	uint64_t *kmer_array_1;
	kmer_array_1=(uint64_t *)malloc(sizeof(uint64_t)*unitnum);
	memset(kmer_array_1,0,sizeof(uint64_t)*unitnum);

	//5.2计算合并次数
	uint32_t N_merge;
	uint32_t N_block;
	N_block=ceil((double)kN/blocksize);
	N_merge=ceil(log(N_block)/log(2));

	uint32_t N_block_tmp=N_block;
	uint64_t *p_start_kmer,*p_des_kmer,*p_kmer_tmp;
	uint32_t *p_start_id,*p_des_id,*p_id_tmp;
	p_start_kmer = kmer_array;
	p_des_kmer = kmer_array_1;
	p_start_id = NULL;
	p_des_id = NULL;
	uint32_t last_kmer;
	for(uint32_t i=0;i<N_merge;i++)
	{
		uint32_t Size_block_loop=blocksize*pow(2,i);
		uint32_t Size_block_loop_p2;

		if(thread_num==1)
		{
			for(uint32_t j=0;j<N_block_tmp/2;j++)
			{

				uint64_t *p1,*p2,*p3;

				p1=p_start_kmer+(2*j)*Size_block_loop*unitperkmer;
				p2=p_start_kmer+(2*j+1)*Size_block_loop*unitperkmer;
				p3=p_des_kmer+(2*j)*Size_block_loop*unitperkmer;


				if(j == (N_block_tmp/2 - 1))
				{
					uint64_t last_p2 = (p2-p_start_kmer)/unitperkmer;
					if(last_p2 + Size_block_loop > kN)
					{
						Size_block_loop_p2 = kN-last_p2;
					}
					else
					{
						Size_block_loop_p2=Size_block_loop;
					}
				}
				else
				{
					Size_block_loop_p2=Size_block_loop;
				}

				uint64_t cnt1 = 0, cnt2 = 0, memcp_n = 0;
				while(cnt1 < Size_block_loop && cnt2 < Size_block_loop_p2)
				{
					if(cmp256BitKmer(p1,p2,unitperkmer)==0)
					{
						cnt1++;
						kmercpy(p3,p1,unitperkmer);
						p1+=unitperkmer;
						p3+=unitperkmer;
					}
					else
					{
						cnt2++;
						kmercpy(p3,p2,unitperkmer);
						p2+=unitperkmer;
						p3+=unitperkmer;
					}
				}
				if(cnt2==Size_block_loop_p2)
				{
					memcp_n = (Size_block_loop - cnt1) ;
					memcpy(p3,p1,memcp_n*unitperkmer * sizeof(uint64_t));
				}
				else
				{
					memcp_n = (Size_block_loop_p2 - cnt2) ;
					memcpy(p3,p2,memcp_n *unitperkmer *  sizeof(uint64_t));
				}
			}
			if(N_block_tmp%2!=0)
			{
				uint32_t offset = Size_block_loop*(N_block_tmp/2)*2;
				last_kmer = kN - offset;
				uint64_t *p_des_kmer_tmp = p_des_kmer + offset*unitperkmer;
				uint64_t *p_start_kmer_tmp = p_start_kmer + offset*unitperkmer;
				memcpy(p_des_kmer_tmp,p_start_kmer_tmp,last_kmer*unitperkmer*sizeof(uint64_t));
			}
		}
		else
		{
			pthread_t* t;
			t=(pthread_t*)malloc(sizeof(pthread_t)*thread_num);
			struct para_merge * p_para;
			p_para=(struct para_merge *)malloc(sizeof(struct para_merge)*thread_num);
			for(uint32_t i_para=0;i_para<thread_num;i_para++)
			{
				p_para[i_para].p_start_kmer=p_start_kmer;
				p_para[i_para].p_des_kmer=p_des_kmer;
				p_para[i_para].p_start_id=p_start_id;
				p_para[i_para].p_des_id=p_des_id;
				p_para[i_para].blocksize=Size_block_loop;
				p_para[i_para].thread_num=thread_num;
				p_para[i_para].thread_id=i_para;
				p_para[i_para].ukN=kN;
				p_para[i_para].unitperkmer=unitperkmer;
				p_para[i_para].N_block=N_block_tmp;
				if(pthread_create(t+i_para, NULL, merge, (void*)(p_para+i_para))!=0)
				{
					cout << "error!" << endl;
				}
			}
			for(uint32_t i_para=0;i_para<thread_num;i_para++)
			{
				pthread_join(t[i_para], NULL);
			}
			free(p_para);
			free(t);
		}

		N_block_tmp=ceil((double)N_block_tmp/2);
		p_kmer_tmp=p_des_kmer;
		p_des_kmer=p_start_kmer;
		p_start_kmer=p_kmer_tmp;

		p_id_tmp=p_des_id;
		p_des_id=p_start_id;
		p_start_id=p_id_tmp;
	}
	cout << "Sort over !" << endl;

	uint64_t *ptr1,*ptr2;
	ptr1 = p_start_kmer;
	ptr2 = p_start_kmer + unitperkmer;
	uint32_t unequ_cnt = 0;

	while(1)
	{
		if(ptr2 == p_start_kmer + kN * unitperkmer)
		{
			break;
		}
		if(cmp256BitKmer(ptr1,ptr2,unitperkmer)==0)
		{
			ptr1 += unitperkmer;
			ptr2 += unitperkmer;
		}
		else
		{
			unequ_cnt++;
//			cout << ptr1 - p_start_kmer << endl;
			ptr1 += unitperkmer;
			ptr2 += unitperkmer;
		}

	}
	cout << "unequ_cnt : " << unequ_cnt << endl;

	char outputFileName[32] = {0};
	sprintf(outputFileName,"%s%s",p_file_path,"_S");
	FILE * fp_output;
	fp_output = fopen(outputFileName,"wb");
	fwrite(p_start_kmer,sizeof(uint64_t)*unitnum,1,fp_output);
	fclose(fp_output);
	cout << "写入完毕" << endl;
//	free(kmer_array_1);
//	free(kmer_array);
//	free(hashvalue_tmp);
}
uint64_t * SortFile_umers_for_split(uint64_t unitnum,uint64_t *kmer_array,uint32_t thread_num,uint8_t unitperkmer){
    uint64_t kN = unitnum/unitperkmer;
    uint64_t * const hashvalue_tmp = (uint64_t *)malloc(sizeof(uint64_t) * unitperkmer);
    cout << "unitnum:" << unitnum << endl;
    cout << "ByteNum:" << unitnum * sizeof(uint64_t) << endl;

    //4)分块排序
    //4.2)分块排序
    uint64_t * cmp_p1,*cmp_p2;
    cmp_p1 = kmer_array;
    cmp_p2 = kmer_array + unitperkmer;
    while(cmp_p1 != kmer_array + kN * unitperkmer && cmp_p2 != kmer_array + kN * unitperkmer)
    {
        if(cmp256BitKmer(cmp_p1,cmp_p2,unitperkmer)==1)
        {
            kmercpy(hashvalue_tmp,cmp_p1,unitperkmer);
            kmercpy(cmp_p1,cmp_p2,unitperkmer);
            kmercpy(cmp_p2,hashvalue_tmp,unitperkmer);
        }
        cmp_p1 += 2*unitperkmer;
        cmp_p2 += 2*unitperkmer;
    }

    //5)整体排序
    //5.1分配轮换空间
    uint32_t blocksize = (uint32_t)pow(2,1);
    uint64_t *kmer_array_1;
    kmer_array_1=(uint64_t *)malloc(sizeof(uint64_t)*unitnum);
    memset(kmer_array_1,0,sizeof(uint64_t)*unitnum);

    //5.2计算合并次数
    uint32_t N_merge;
    uint32_t N_block;
    N_block=ceil((double)kN/blocksize);
    N_merge=ceil(log(N_block)/log(2));

    uint32_t N_block_tmp=N_block;
    uint64_t *p_start_kmer,*p_des_kmer,*p_kmer_tmp;
    uint32_t *p_start_id,*p_des_id,*p_id_tmp;
    p_start_kmer = kmer_array;
    p_des_kmer = kmer_array_1;
    p_start_id = NULL;
    p_des_id = NULL;
    uint32_t last_kmer;
    for(uint32_t i=0;i<N_merge;i++)
    {
        uint32_t Size_block_loop=blocksize*pow(2,i);
        uint32_t Size_block_loop_p2;

        if(thread_num==1)
        {
            for(uint32_t j=0;j<N_block_tmp/2;j++)
            {

                uint64_t *p1,*p2,*p3;

                p1=p_start_kmer+(2*j)*Size_block_loop*unitperkmer;
                p2=p_start_kmer+(2*j+1)*Size_block_loop*unitperkmer;
                p3=p_des_kmer+(2*j)*Size_block_loop*unitperkmer;


                if(j == (N_block_tmp/2 - 1))
                {
                    uint64_t last_p2 = (p2-p_start_kmer)/unitperkmer;
                    if(last_p2 + Size_block_loop > kN)
                    {
                        Size_block_loop_p2 = kN-last_p2;
                    }
                    else
                    {
                        Size_block_loop_p2=Size_block_loop;
                    }
                }
                else
                {
                    Size_block_loop_p2=Size_block_loop;
                }

                uint64_t cnt1 = 0, cnt2 = 0, memcp_n = 0;
                while(cnt1 < Size_block_loop && cnt2 < Size_block_loop_p2)
                {
                    if(cmp256BitKmer(p1,p2,unitperkmer)==0)
                    {
                        cnt1++;
                        kmercpy(p3,p1,unitperkmer);
                        p1+=unitperkmer;
                        p3+=unitperkmer;
                    }
                    else
                    {
                        cnt2++;
                        kmercpy(p3,p2,unitperkmer);
                        p2+=unitperkmer;
                        p3+=unitperkmer;
                    }
                }
                if(cnt2==Size_block_loop_p2)
                {
                    memcp_n = (Size_block_loop - cnt1) ;
                    memcpy(p3,p1,memcp_n*unitperkmer * sizeof(uint64_t));
                }
                else
                {
                    memcp_n = (Size_block_loop_p2 - cnt2) ;
                    memcpy(p3,p2,memcp_n *unitperkmer *  sizeof(uint64_t));
                }
            }
            if(N_block_tmp%2!=0)
            {
                uint32_t offset = Size_block_loop*(N_block_tmp/2)*2;
                last_kmer = kN - offset;
                uint64_t *p_des_kmer_tmp = p_des_kmer + offset*unitperkmer;
                uint64_t *p_start_kmer_tmp = p_start_kmer + offset*unitperkmer;
                memcpy(p_des_kmer_tmp,p_start_kmer_tmp,last_kmer*unitperkmer*sizeof(uint64_t));
            }
        }
        else
        {
            pthread_t* t;
            t=(pthread_t*)malloc(sizeof(pthread_t)*thread_num);
            struct para_merge * p_para;
            p_para=(struct para_merge *)malloc(sizeof(struct para_merge)*thread_num);
            for(uint32_t i_para=0;i_para<thread_num;i_para++)
            {
                p_para[i_para].p_start_kmer=p_start_kmer;
                p_para[i_para].p_des_kmer=p_des_kmer;
                p_para[i_para].p_start_id=p_start_id;
                p_para[i_para].p_des_id=p_des_id;
                p_para[i_para].blocksize=Size_block_loop;
                p_para[i_para].thread_num=thread_num;
                p_para[i_para].thread_id=i_para;
                p_para[i_para].ukN=kN;
                p_para[i_para].unitperkmer=unitperkmer;
                p_para[i_para].N_block=N_block_tmp;
                if(pthread_create(t+i_para, NULL, merge, (void*)(p_para+i_para))!=0)
                {
                    cout << "error!" << endl;
                }
            }
            for(uint32_t i_para=0;i_para<thread_num;i_para++)
            {
                pthread_join(t[i_para], NULL);
            }
            free(p_para);
            free(t);
        }

        N_block_tmp=ceil((double)N_block_tmp/2);
        p_kmer_tmp=p_des_kmer;
        p_des_kmer=p_start_kmer;
        p_start_kmer=p_kmer_tmp;

        p_id_tmp=p_des_id;
        p_des_id=p_start_id;
        p_start_id=p_id_tmp;
    }
    cout << "Sort over !" << endl;

    uint64_t *ptr1,*ptr2;
    ptr1 = p_start_kmer;
    ptr2 = p_start_kmer + unitperkmer;
    uint32_t unequ_cnt = 0;

    while(1)
    {
        if(ptr2 == p_start_kmer + kN * unitperkmer)
        {
            break;
        }
        if(cmp256BitKmer(ptr1,ptr2,unitperkmer)==0)
        {
            ptr1 += unitperkmer;
            ptr2 += unitperkmer;
        }
        else
        {
            unequ_cnt++;
//            cout << ptr1 - p_start_kmer << endl;
            ptr1 += unitperkmer;
            ptr2 += unitperkmer;
        }

    }
    cout << "unequ_cnt : " << unequ_cnt << endl;

    return p_start_kmer;
    free(kmer_array_1);
    free(kmer_array);
    free(hashvalue_tmp);
}