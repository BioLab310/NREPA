//
// Created by lin on 2021/11/18.
//
#include "FM-index.h"
void fm_main(int argc, char **argv){
    struct build_para para;
    char c;
    uint32_t k_len;
    for(int32_t i=2;i<argc;i=i+2)
    {
        if(argv[i][0]=='-'&&argv[i][1]=='h')	//path of reference
        {
            para.ref_path=argv[i+1];
        }
        else if(argv[i][0]=='-'&&argv[i][1]=='s')//sa_gap
        {
            para.sa_gap=atoi(argv[i+1]);
        }
        else if(argv[i][0]=='-'&&argv[i][1]=='o')//sa_gap
        {
            para.occ_gap=atoi(argv[i+1]);
        }
        else if(argv[i][0]=='-'&&argv[i][1]=='l')//level
        {
            para.level=atoi(argv[i+1]);
        }
        else if(argv[i][0]=='-'&&argv[i][1]=='k')//k_mer
        {
            k_len=atoi(argv[i+1]);
        }
        else if(argv[i][0]=='-'&&argv[i][1]=='t')//thread_num
        {
            para.thread_num=atoi(argv[i+1]);
        }
    }

    char *dir = "./index";
    if(access(dir,0) == -1)
    {
        printf("%s not exist.\n",dir);
        mkdir(dir, 0755);
    }
    build(para);
//    build_occA(para);
}

void build(struct build_para &para_build){
    char * ref_path=para_build.ref_path;
    uint32_t thread_log=para_build.thread_num;
    uint32_t sa_gap=para_build.sa_gap;
    uint32_t level=para_build.level;
    uint64_t max_len=para_build.max_len;

    FILE * out_C = NULL;
    FILE * out_B = NULL;
    out_C=fopen("./index/C","wb+");
    out_B=fopen("./index/B","wb+");

    char * ref;
    uint32_t len;

    ReadSeq(&ref,&len,ref_path);

    struct para_cmp para;
    para.ref=ref;
    para.ref_len=len;
    para.max_len=max_len;

    char *b;
    b=(char*)malloc(sizeof(char)*(len+1));

    uint32_t C[5]={0};
    for(uint32_t i=0;i<len;i++)
    {
        switch(ref[i])
        {
            case 'A':
                C[0]++;
                break;
            case 'C':
                C[1]++;
                break;
            case 'G':
                C[2]++;
                break;
            case 'T':
                C[3]++;
                break;
            case '$':
                C[4]++;
                break;
            default:
                break;
        }
    }
    cout << "$的数目：" << C[4] << endl;
    fwrite(C,sizeof(uint32_t),5,out_C);

    uint32_t dollar_pos;
    b[0]=ref[len-1];
    uint32_t b_len_cur=1;

    uint32_t block_num=pow(8,level);
    char * flag;
    flag=(char*)malloc(sizeof(char)*level);
    uint32_t block_start=1;

    for(uint32_t l=0;l<block_num;l++)
    {
        uint32_t plate=7;
        uint32_t sign = 0;
        for(uint32_t k=0;k<level;k++)
        {
            uint32_t tmp=l>>(3*(level-1-k));
            switch (tmp&plate)
            {
                case 0:
                    flag[k]='$';
                    break;
                case 1:
                    flag[k]='A';
                    break;
                case 2:
                    flag[k]='C';
                    break;
                case 3:
                    flag[k]='G';
                    break;
                case 4:
                    flag[k]='T';
                    break;
                default:
                    sign = 1;
                    break;
            }
            if(sign == 1)
                break;
        }
        if(sign == 1)
            continue;
        if(flag[0]=='$' && flag[1]=='$')
            continue;
        cout << flag << endl;

        struct timeval tvs,tve;
        gettimeofday(&tvs,NULL);

        struct Node * root;
        root=NodeInitial();

        for(uint32_t i=0;i<len;i++)
        {
            uint32_t f=1;
            for(uint32_t x=0;x<level;x++)
            {
                if(i+x<len)
                {
                    if(ref[i+x]!=flag[x])
                    {
                        f=0;
                        break;
                    }
                }
                else
                {
                    if(flag[x]!='A')
                    {
                        f=0;
                        break;
                    }
                }
            }
            if(f==1)
            {
                char tmp;
                if(i==0)
                {
                    tmp='#';
                }
                else
                {
                    tmp=para.ref[i-1];
                }
                insert(&root,root, i,tmp,para);
            }
        }

        gettimeofday(&tve,NULL);
        double span = tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;
        cout <<"time of construct B+tree is: "<<span<<endl;

        struct Node * mini;
        mini=Find_Minimal_Node(root);

        while(mini!=NULL)
        {
            for(uint8_t i=0;i<mini->nodesize;i++)
            {
                b[b_len_cur]=mini->b[i];
                if(mini->b[i]=='#')
                {
                    dollar_pos=b_len_cur;
                }

                b_len_cur++;
            }
            mini=mini->brother;
        }
        treeDestroy(root);
    }
    fwrite(&para_build.sa_gap,sizeof(uint32_t),1,out_C);
    fwrite(&para_build.occ_gap,sizeof(uint32_t),1,out_C);
    fwrite(b,sizeof(char),len+1,out_B);

    free(b);
    fclose(out_C);
    fclose(out_B);
}

void build_occA(struct build_para &para_build)
{
    char * b = NULL;
    uint32_t len = 0;

    FILE * B_file = NULL;
    B_file = fopen("./index/B","r");

    FILE * C_file = NULL;
    C_file = fopen("./index/C","r");

    FILE * out_OCCA = NULL;
    out_OCCA = fopen("./index/OCCA","wb+");

    uint64_t C[5]={0};
    if(C_file){
        fread(C,sizeof(uint64_t),5,C_file);
    }

    if(B_file)//打开文件一定要判断是否成功
    {
        fseek(B_file,0,SEEK_END);//将文件内部的指针指向文件末尾
        len = ftell(B_file);//获取文件长度，（得到文件位置指针当前位置相对于文件首的偏移字节数）
        rewind(B_file);//将文件内部的指针重新指向一个流的开头
        b = (char *) malloc(len*sizeof(char));//申请内存空间，len*sizeof(char)是为了更严谨，16位上char占一个字符，其他机器上可能变化

        //用malloc申请的内存是没有初始值的，如果不赋值会导致写入的时候找不到结束标志符而出现内存比实际申请值大，写入数据后面跟随乱码的情况
        memset(b,0,len*sizeof(char));//将内存空间都赋值为‘\0’

        fread(b,sizeof(char),len,B_file);
    }
    cout << "$ num is " << C[4] << endl;
    uint32_t count = 0;
    uint32_t *a;
    a = (uint32_t*)malloc(sizeof(uint32_t)*4);
    memset(a,0,sizeof(uint32_t)*4);
    for(uint32_t i = 0; i <= len; i++)
    {
        if(i != 0)
        {
            switch(b[i-1])
            {
                case 'A':
                    a[0]++;
                    break;
                case 'C':
                    a[1]++;
                    break;
                case 'G':
                    a[2]++;
                    break;
                case 'T':
                    a[3]++;
                    break;
                default:
                    break;
            }
        }
        if(i >= C[4] && ((i - C[4]) % para_build.occ_gap) == 0)
        {
//			fprintf(out_OCCA, "%b ", a[0]);
//			fprintf(out_OCCA, "%b ", a[1]);
//			fprintf(out_OCCA, "%b ", a[2]);
//			fprintf(out_OCCA, "%b ", a[3]);
//            cout << a[0] << endl;
//            cout << a[1] << endl;
//            cout << a[2] << endl;
//            cout << a[3] << endl;
            fwrite(a,sizeof(uint32_t)*4,1,out_OCCA);
            count++;
        }
    }
    cout << "occa length is " << count << endl;
    free(a);
    free(b);
    fclose(out_OCCA);
    fclose(B_file);
}
