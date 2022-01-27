/*
 * inputFMindex.cpp
 *
 *  Created on: Apr 15, 2021
 *      Author: bio
 */

#include "inputFMindex.h"

void *read_binfile(char *path)
{
	void *ptr;
	FILE * file = NULL;
	file = fopen(path,"r");
	uint32_t len;
	uint32_t pathlen;
	pathlen = strlen(path);
	if(path[pathlen-1] == 'B')
	{
		if(file)
		{
			fseek(file,0,SEEK_END);//将文件内部的指针指向文件末尾
			len = ftell(file);//获取文件长度，（得到文件位置指针当前位置相对于文件首的偏移字节数）
			rewind(file);//将文件内部的指针重新指向一个流的开头
			ptr = malloc(len*sizeof(char));//申请内存空间，len*sizeof(char)是为了更严谨，16位上char占一个字符，其他机器上可能变化

			//用malloc申请的内存是没有初始值的，如果不赋值会导致写入的时候找不到结束标志符而出现内存比实际申请值大，写入数据后面跟随乱码的情况
			memset(ptr,0,len*sizeof(char));//将内存空间都赋值为‘\0’

			fread(ptr,sizeof(char),len,file);
		}
	}
	else
	{
		if(file)
		{
			fseek(file,0,SEEK_END);//将文件内部的指针指向文件末尾
			len = ftell(file);//获取文件长度，（得到文件位置指针当前位置相对于文件首的偏移字节数）
			rewind(file);//将文件内部的指针重新指向一个流的开头
			ptr = malloc(len);//申请内存空间，len*sizeof(char)是为了更严谨，16位上char占一个字符，其他机器上可能变化
			//用malloc申请的内存是没有初始值的，如果不赋值会导致写入的时候找不到结束标志符而出现内存比实际申请值大，写入数据后面跟随乱码的情况
			memset(ptr,0,len);//将内存空间都赋值为‘\0’

			fread(ptr,1,len,file);
		}
	}
	fclose(file);
	return (void *)(ptr);
}

void read_bfile2index(const char *path, sFMindex &FMidx, uint32_t num)
{
	//num is the number of file :  SA1 SA2 SA3...
	FILE * file = NULL;
	char filename[64] = {0};
	size_t len;
	if(num == 0)
	{
		sprintf(filename,"%s%s",path,"/B");
	}
	else
	{
		sprintf(filename,"%s%s%d",path,"/B",num);
	}
	file = fopen(filename,"r");
	if(file)
	{
		fseek(file,0,SEEK_END);//将文件内部的指针指向文件末尾
		len = ftell(file);//获取文件长度，（得到文件位置指针当前位置相对于文件首的偏移字节数）
		rewind(file);//将文件内部的指针重新指向一个流的开头
		FMidx.b = (char *)malloc(len+1*sizeof(char));//申请内存空间，len*sizeof(char)是为了更严谨，16位上char占一个字符，其他机器上可能变化
		//用malloc申请的内存是没有初始值的，如果不赋值会导致写入的时候找不到结束标志符而出现内存比实际申请值大，写入数据后面跟随乱码的情况
		memset(FMidx.b,0,len+1);//将内存空间都赋值为‘\0’

		fread(FMidx.b, len, 1, file);
		fclose(file);
	}
	if(num == 0)
	{
		sprintf(filename,"%s%s",path,"/C");
	}
	else
	{
		sprintf(filename,"%s%s%d",path,"/C",num);
	}
	file = fopen(filename,"r");
	if(file)
	{
		//C file include 6 uint32_t type date respective represent number of A,C,G,T and sa_gap occ_gap
		uint32_t totallen = sizeof(uint32_t) * 5;
		FMidx.c = (uint32_t *)malloc(totallen);
		//用malloc申请的内存是没有初始值的，如果不赋值会导致写入的时候找不到结束标志符而出现内存比实际申请值大，写入数据后面跟随乱码的情况
		memset(FMidx.c,0,totallen);//将内存空间都赋值为‘\0’
		size_t readnum = fread(FMidx.c + 1,sizeof(uint32_t),4,file);
		if(readnum != 4)
		{
			exit(4);
		}
		for(int i = 2; i < 5; ++i)
		{
			FMidx.c[i] += FMidx.c[i-1];
		}
		fread(&FMidx.sa_gap,sizeof(FMidx.sa_gap),1,file);
		fread(&FMidx.occ_gap,sizeof(FMidx.occ_gap),1,file);
		fclose(file);
	}
	if(num == 0)
	{
		sprintf(filename,"%s%s",path,"/OCCA");
	}
	else
	{
		sprintf(filename,"%s%s%d",path,"/OCCA",num);
	}
	file = fopen(filename,"r");
	if(file)
	{
		fseek(file,0,SEEK_END);//将文件内部的指针指向文件末尾
		len = ftell(file);//获取文件长度，（得到文件位置指针当前位置相对于文件首的偏移字节数）
		rewind(file);//将文件内部的指针重新指向一个流的开头
		FMidx.occa = (uint32_t *)malloc(len);//申请内存空间，len*sizeof(char)是为了更严谨，16位上char占一个字符，其他机器上可能变化
		//用malloc申请的内存是没有初始值的，如果不赋值会导致写入的时候找不到结束标志符而出现内存比实际申请值大，写入数据后面跟随乱码的情况
		memset(FMidx.occa,0,len);//将内存空间都赋值为‘\0’

		fread(FMidx.occa,len,1,file);
		fclose(file);
		FMidx.occ_num = len / 4;
	}
	if(num == 0)
	{
		sprintf(filename,"%s%s",path,"/SA");
	}
	else
	{
		sprintf(filename,"%s%s%d",path,"/SA",num);
	}
	file = fopen(filename,"r");
	if(file)
	{
		fseek(file,0,SEEK_END);//将文件内部的指针指向文件末尾
		len = ftell(file);//获取文件长度，（得到文件位置指针当前位置相对于文件首的偏移字节数）
		rewind(file);//将文件内部的指针重新指向一个流的开头
		FMidx.sa = (uint32_t *)malloc(len);//申请内存空间，len*sizeof(char)是为了更严谨，16位上char占一个字符，其他机器上可能变化
		//用malloc申请的内存是没有初始值的，如果不赋值会导致写入的时候找不到结束标志符而出现内存比实际申请值大，写入数据后面跟随乱码的情况
		memset(FMidx.sa,0,len);//将内存空间都赋值为‘\0’

		fread(FMidx.sa,len,1,file);
		fclose(file);
		FMidx.sa_num = len / 4;
	}
}

void free_FMindex(sFMindex &FMidx)
{
	if(FMidx.b)
	{
		free(FMidx.b);
		FMidx.b = NULL;
	}
	if(FMidx.c)
	{
		free(FMidx.c);
		FMidx.c = NULL;
	}
	if(FMidx.occa)
	{
		free(FMidx.occa);
		FMidx.occa = NULL;
	}
	if(FMidx.sa)
	{
		free(FMidx.sa);
		FMidx.sa = NULL;
	}
}




