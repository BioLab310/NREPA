//
// Created by lin on 21-10-25.
//
#include "utils.h"

void utils_main(int argc, char **argv){

    cout << "start merge" << endl;

    char **p_ref_path; //ref參考基因 用二维指针 因为可能读入多个基因组
    //基因组的个数和基因组名字长度 暂时定为20
    char *sort_file_path;//要排序的文件
    uint32_t thread_num;//线程数目
    uint32_t k_mer;
    p_ref_path = (char**) malloc(sizeof(char*) * 100);
    *p_ref_path = (char*) malloc(sizeof(char) * 100);
    uint32_t f_count = 0; //读取的基因组的个数

    for (int32_t i = 2; i < argc; i++) {
        if (argv[i][0] == '-' && argv[i][1] == 'R') {
            p_ref_path[f_count++] = argv[i + 1];
            for (int32_t j = 2; j <= argc - 3; j++) {
                if (argv[i + j][0] != '-') {
                    p_ref_path[f_count++] = argv[i + j];
                } else {
                    break;
                }
            }
            merge(f_count,p_ref_path);
        }
        if (argv[i][0] == '-' && argv[i][1] == 'c') {
//            uint64_t copy_size = 2898295643;
//            copy_utils(copy_size);
            writeArray();
//            readArray("a");
//            break;
        }
        if (argv[i][0] == '-' && argv[i][1] == 't') {
            thread_num = atoi(argv[i + 1]);
        }
        if (argv[i][0] == '-' && argv[i][1] == 'u') {
            sort_file_path = argv[i + 1];
            SortFile_umers_forB(sort_file_path, thread_num, k_code_size);
            break;
        }
        if (argv[i][0] == '-' && argv[i][1] == 'p') {
            split_detail(argc,argv);
            break;
        }
        if (argv[i][0] == '-' && argv[i][1] == '$') {;
            char *path = argv[i + 1];
            cout << "计算"<< path << "中$数目" << endl;
            splitNum_v2(path);
            break;
        }
    }
    cout << "end merge" << endl;
}

void split_detail(int argc, char **argv){
    char *ref_path;
    char *sav_path;
    uint64_t start;
    uint64_t span;

    for (int32_t i = 3; i < argc; i++){
        if(argv[i][0] == '-' && argv[i][1] == 'r'){
            ref_path = argv[i + 1];
        }
        if(argv[i][0] == '-' && argv[i][1] == 's'){
            start = atoi(argv[i + 1]);
        }
        if(argv[i][0] == '-' && argv[i][1] == 'n'){
            span = atoi(argv[i + 1]);
        }
        if(argv[i][0] == '-' && argv[i][1] == 'v'){
            sav_path = argv[i + 1];
        }
    }
    split_compress_ref(ref_path,start,span,sav_path);
}

void copy_utils(uint64_t copy_size){
    cout << "utils for copy correct char" << endl;
    char *p_ref;
    uint_fast64_t ref_length;
    char *p_ref_path = "compress_ref.fq";
    ReadSeq(&p_ref, &ref_length, p_ref_path);
    char *copy;
    copy = (char *)malloc(sizeof(char)*copy_size);
    for(uint64_t i=0;i<copy_size;i++){
        copy[i]=p_ref[i];
    }

    //打开之前保存的压缩文件
    char *p_ref_compress;
    uint64_t p_ref_compress_length;
    //写入compress_ref
    FILE *write_ref_file = fopen("compress_ref_v2.fq", "w");
    if(!write_ref_file){
        cout << "can't find compress_ref.fq" << endl;
    }
    fwrite(copy, sizeof(char), copy_size, write_ref_file);
    fclose(write_ref_file);
    cout << "end copy" << endl;

    cout << "test my copy" << endl;
    FILE *fin = fopen("compress_ref_v2.fq","r");
    if(!fin){
        cout << "can't find compress_ref_v2.fq" << endl;
    }
    fseek(fin, 0, 2);
    //计算序列总长度
    p_ref_compress_length = ftell(fin);
    //指针定位到文件头
    fseek(fin, 0, 0);
    p_ref_compress = (char *)malloc(sizeof(char)*p_ref_compress_length);
    fread(p_ref_compress, sizeof(char), p_ref_compress_length, fin);
    fclose(fin);
    cout << "文件保存的compress length：" << p_ref_compress_length << endl;
    cout << "last char:" << p_ref_compress[p_ref_compress_length-1] << endl;
}

void merge(uint32_t f_count,char **p_ref_path){
    //输出多个基因组合并结果 MAXSIZE是合并后的基因的长度 不同服务器它的长度不一样
    uint64_t maxsize = MAXSIZE;
    char *mul_ref;
    mul_ref = (char*) malloc(sizeof(char) * maxsize);
    uint64_t mul_ref_length;

    //每一条基因和其长度
    char *p_ref;
    uint64_t ref_length;
    for (uint32_t i = 0; i < f_count; i++) {
        ReadSeq(&p_ref, &ref_length, p_ref_path[i]);
        //将多个基因组拼接在一起，最后输出合并后的基因和长度
        strcat(mul_ref, p_ref);
    }
    mul_ref_length = strlen(mul_ref);

    //写入compress_ref
    FILE *write_ref_file = fopen("merge.fq", "w");
    if(!write_ref_file){
        cout << "can't find merge.fq" << endl;
    }
    fwrite(mul_ref, sizeof(char), mul_ref_length, write_ref_file);
    fclose(write_ref_file);
    cout << "write length:" << mul_ref_length << endl;

    cout << "check out write" << endl;
    char *check_ref;
    uint64_t check_length;
    FILE *fin = fopen("merge.fq","r");
    if(!fin){
        cout << "can't find compress_ref.fq" << endl;
    }
    fseek(fin, 0, 2);
    //计算序列总长度
    check_length = ftell(fin);
    //指针定位到文件头
    fseek(fin, 0, 0);
    check_ref = (char *)malloc(sizeof(char)*check_length);
    fread(check_ref, sizeof(char), check_length, fin);
    fclose(fin);
    cout << "read length:" << check_length << endl;
}

void split_compress_ref(char *p_ref_path,uint64_t start,uint64_t span,char *p_sav_path){
    cout << "start" << endl;
    char *p_ref;
    uint64_t ref_length;
    uint64_t last_len;
    start = 0;
    ReadSeq(&p_ref, &ref_length, p_ref_path);
    char *split;
    split = (char *)malloc(sizeof(char)*ref_length);
    FILE *write_ref_file = fopen(p_sav_path, "w");
    if(!write_ref_file){
        cout << "can't find compress_ref_v2.fq" << endl;
    }
    for(uint64_t i=start;i<ref_length;i++){
        split[i-start] = p_ref[i];
        if(i>start + (span*100000000) && p_ref[i]=='$'){
            split[i-start]=p_ref[i];
            last_len = i - start + 1;
            fwrite(split, sizeof(char), last_len, write_ref_file);
            fclose(write_ref_file);
            break;
        }
        if(i==ref_length-1){
            split[i-start]=p_ref[i];
            last_len = i - start + 1;
            fwrite(split, sizeof(char), last_len, write_ref_file);
            fclose(write_ref_file);
            break;
        }
    }
    cout << "数组最后一个序号是" << start + last_len - 1 << endl;
    cout << "最后一个字符是" << split[last_len - 1] << endl;
}
void splitNum_v2(char *p_ref_path){

    //每个子串压缩后的p_ref
    char *ref;
    uint32_t ref_length = 0;
    ReadSeq(&ref, &ref_length, p_ref_path);

    uint32_t C[5]={0};

    for(uint32_t i = 0; i < ref_length; i++){

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
    cout << C[0] << endl;
    cout << C[1] << endl;
    cout << C[2] << endl;
    cout << C[3] << endl;
    cout << C[4] << endl;
}

void writeArray(){
    FILE *B_file = fopen("B", "r");
    uint64_t len;
    char *b = NULL;
    if(!B_file) {
        cout << "have not B" << endl;
    }
    fseek(B_file,0,SEEK_END);//将文件内部的指针指向文件末尾
    len = ftell(B_file);//获取文件长度，（得到文件位置指针当前位置相对于文件首的偏移字节数）
    rewind(B_file);//将文件内部的指针重新指向一个流的开头
    b = (char *) malloc(len*sizeof(char));//申请内存空间，len*sizeof(char)是为了更严谨，16位上char占一个字符，其他机器上可能变化
    //用malloc申请的内存是没有初始值的，如果不赋值会导致写入的时候找不到结束标志符而出现内存比实际申请值大，写入数据后面跟随乱码的情况
    memset(b,0,len*sizeof(char));//将内存空间都赋值为‘\0’
    fread(b,sizeof(char),len,B_file);
    cout << len << endl;

    uint32_t C[5]={0};
    FILE * C_file = NULL;
    C_file = fopen("C","r");

    if(C_file){
        fread(C,sizeof(uint32_t),5,C_file);
    }
    cout << C[0] << endl;
    cout << C[1] << endl;
    cout << C[2] << endl;
    cout << C[3] << endl;
    cout << C[4] << endl;

    uint32_t *a;
    a = (uint32_t*)malloc(sizeof(uint32_t)*4);
    memset(a,0,sizeof(uint32_t)*4);
    uint32_t *sav = (uint32_t *)malloc(sizeof(uint32_t)*(len-C[4]+1)*4);
    uint64_t j_count = 0;
    cout << "开始计算" << endl;
    for(uint32_t i = 0; i <= len; i++){
//        if(i%100000==0){
//            cout << i << endl;
//        }
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

        if(i >= C[4] && ((i - C[4]) % 1) == 0)
        {
            for(uint32_t j = 0;j<4;j++){
                sav[j_count++] = a[j];
            }

        }

    }
    FILE *save_binary_file = fopen("OCCA", "wb");
    fwrite(sav, sizeof(uint32_t), (len-C[4]+1)*4, save_binary_file);
    fclose(B_file);
    fclose(C_file);
    fclose(save_binary_file);
}

void readArray(char *outputFileName){
    FILE *read_binary_file = fopen("B", "r");
    //指针定位到文件尾
    fseek(read_binary_file, 0, 2);
    //计算一共几个64位
    uint64_t total_size = ftell(read_binary_file)/sizeof(uint64_t);
    cout << total_size << endl;
    //指针定位到文件头
    fseek(read_binary_file, 0, 0);

    uint64_t *myArray= (uint64_t*) malloc(
            sizeof(uint64_t) * total_size);
    fread(myArray, sizeof(uint32_t), total_size, read_binary_file);


    fclose(read_binary_file);
//    cout << myArray[3000000000] << endl;
    cout << myArray[total_size-1] << endl;
}