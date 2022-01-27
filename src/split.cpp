//
// Created by lin on 21-10-25.
//
#include "split.h"
void split_main(int argc, char **argv){
    //读入的基因组
    char *p_ref_path;
    //k_mer长度
    uint32_t k_mer;

    for (int32_t i = 2; i < argc; i++){
        if (argv[i][0] == '-' && argv[i][1] == 'R') {
            p_ref_path = argv[i + 1];
        }
        if (argv[i][0] == '-' && argv[i][1] == 'k') {
            k_mer = atoi(argv[i + 1]);
            //k_mer长度决定保存每个k_mer的文件大小
            k_code_size = kmer_codeSize(k_mer);
        }
    }
    splitRef(p_ref_path,k_mer);
}

void splitRef(char *p_split_path,uint32_t k_mer){
    cout << "splitRef start" << endl;
    map<uint64_t,uint64_t> my_count;
    cout << "计算"<< p_split_path << "中$数目" << endl;
    my_count = splitNum(p_split_path);

    vector<string> res;
    vector<string> res_v;
    //数组的数目
    cout << "得到数组的数目" << endl;
    uint64_t count = res_split(res,my_count,p_split_path);
    cout << "数组的数目是：" << count << endl;
    uint64_t *start_kmer;
    uint64_t *end_kmer;
    string s_res = split_merge(res,k_mer,count,&start_kmer,&end_kmer);
    cout << "合并的字符串长度" << s_res.length() << endl;
    ofstream w_file("res_str.fq");
    w_file << s_res;
    w_file.close();
    cout << "splitRef end" << endl;
}

uint64_t res_split(vector<string> &res,map<uint64_t,uint64_t> my_count,char *path){
    string temp;
    ifstream readFile(path);
    getline(readFile,temp);
    uint64_t start=0;
    uint64_t end=0;
    uint64_t count = 0;
    for(uint64_t i=0;i<my_count.size()-1;i++){
        if(my_count.at(i)!=my_count.at(i+1)){
            if(temp[my_count.at(i)]=='$'){
                start=my_count.at(i)+1;
            }else{
                start=my_count.at(i);
            }
            if(temp[my_count.at(i+1)]=='$'){
                end=my_count.at(i+1)-1;
            }else{
                end=my_count.at(i+1);
            }
            if(start!=end){
                res.push_back(temp.substr(start,end-start+1).c_str());
                count++;
            }
        }
    }
    cout << "共保存了" << res.size() << "个数组" << endl;
    return count;
}

string split_merge(vector<string> &res,uint32_t k_mer,uint64_t count,uint64_t **start_kmer,uint64_t **end_kmer){
    cout << "开始合并成一个字符串" << endl;
    string res_str;
    struct bit256KmerPara para_tmp;
    unordered_set<uint64_t> set_index;
    map<uint64_t,uint64_t> map_index;
    vector<string> res_v;

    *start_kmer = (uint64_t *) malloc(sizeof(uint64_t)*count*(k_code_size+1));
    *end_kmer = (uint64_t *) malloc(sizeof(uint64_t)*count*(k_code_size+1));

    uint64_t *start_k_code;
    //定义为3,因为k_mer长度不超过96
    start_k_code = (uint64_t*) malloc(sizeof(uint64_t) * k_code_size);

    uint64_t *end_k_code;
    //定义为3,因为k_mer长度不超过96
    end_k_code = (uint64_t*) malloc(sizeof(uint64_t) * k_code_size);

    /**
 * 开始kmer的code和结束kmer的code以及在数组的位置
 */
    for(uint64_t i=0;i < count;i++){
        memset(start_k_code,0,sizeof(uint64_t) * k_code_size);
        memset(end_k_code,0,sizeof(uint64_t) * k_code_size);
        get_para(&para_tmp, k_mer);
        //分別求开始和结束的k_mer的哈希值
        cal_hash_value_directly_256bit_const(res[i].c_str(), start_k_code, para_tmp);
        cal_hash_value_directly_256bit_const(res[i].c_str()+ res[i].length()-k_mer, end_k_code, para_tmp);
        for(uint32_t k=0;k<k_code_size;k++){
            (*start_kmer)[i*(k_code_size+1)+k] = start_k_code[k];
            (*end_kmer)[i*(k_code_size+1)+k] = end_k_code[k];
        }
        //存储code以及在数组中的位置，k_code_size=3的话应该用4个uint64 第4个保存位置
        (*start_kmer)[i*(k_code_size+1)+k_code_size] = i;
        (*end_kmer)[i*(k_code_size+1)+k_code_size] = i;
    }

    cout << "开始排序" << endl;
    uint64_t * sort_start_kmer = SortFile_umers_for_split(count*(k_code_size+1),*start_kmer,1,k_code_size+1);
    uint64_t * sort_end_kmer = SortFile_umers_for_split(count*(k_code_size+1),*end_kmer,1,k_code_size+1);
    cout << "排序完毕" << endl;

    //计算需要拼接的数组
    uint64_t i=0;
    uint64_t j=0;
    while(i<count&& j<count){
        if(cmp256BitKmer(sort_end_kmer+i*(k_code_size+1),sort_start_kmer+j*(k_code_size+1),k_code_size)==2
           &&sort_end_kmer[i*(k_code_size+1)+k_code_size] != sort_start_kmer[j*(k_code_size+1)+k_code_size]){

            map_index.insert(pair<uint64_t,uint64_t>(sort_end_kmer[i*(k_code_size+1)+k_code_size],sort_start_kmer[j*(k_code_size+1)+k_code_size]));
            set_index.insert(sort_end_kmer[i*(k_code_size+1)+k_code_size]);
            set_index.insert(sort_start_kmer[j*(k_code_size+1)+k_code_size]);
            i++;
            j++;
        } else if(cmp256BitKmer(sort_end_kmer+i*(k_code_size+1),sort_start_kmer+j*(k_code_size+1),k_code_size)==1){
            j++;
        } else if(cmp256BitKmer(sort_end_kmer+i*(k_code_size+1),sort_start_kmer+j*(k_code_size+1),k_code_size)==0){
            i++;
        }else{
            if(cmp256BitKmer(sort_end_kmer+i*(k_code_size+1),sort_start_kmer+(j+1)*(k_code_size+1),k_code_size)==2){
                j++;
            }else{
                i++;
            }
        }
    }

//    cout << "共有" << count - set_index.size() << "个数组无需拼接" << endl;

    cout << "组合无需拼接的数组" << endl;
    for(uint64_t i = 0;i<count;i++){
        if(!set_index.count(i)){
            res_str += res[i];
            res_str += '$';
        }
    }
    cout << "组合完毕" << endl;

    cout << "无需拼接的字符串长度" << res_str.size() << endl;

    cout << "开始拼接数组" << endl;
//    cout << "需要拼接的数组:" << set_index.size() << endl;

    //直接组合map_index里的字符串
    unordered_set<uint64_t> un_need_index;
    vector<uint64_t> keys;
    vector<uint64_t> values;
    for(auto temp1:map_index){
        keys.push_back(temp1.first);
    }
    sort(keys.begin(),keys.end());
    for(auto temp2:map_index){
        values.push_back(temp2.second);
    }
    sort(values.begin(),values.end());

    uint64_t s = 0;
    uint64_t e = 0;
    while (s < keys.size() && e < values.size()){
        if(keys[s] == values[e]){
            un_need_index.insert(values[e]);
            s++;
            e++;
        }else if(keys[s] < values[e]){
            s++;
        } else{
            e++;
        }
    }

    //直接组合map_index里的字符串
    uint64_t merge_count = 0;
    for(auto temp:map_index){
        string end_kmer;
        string start_kmer;
        string merge;
        if(un_need_index.count(temp.first)){
            res_str += res[temp.second];

        }else{
            end_kmer = res[temp.first];
            start_kmer = res[temp.second];
            merge = end_kmer + start_kmer.substr(k_mer);
            res_str += merge;
            merge_count ++;
        }
        res_str += '$';
    }
    cout << "拼接了" << merge_count << "次" << endl;
    cout << "拼接完毕" << endl;
    cout << "组合完毕" << endl;
    return res_str;
}

map<uint64_t,uint64_t> splitNum(char *p_ref_path){

    //每个子串压缩后的p_ref
    FILE *read_ref_file = fopen(p_ref_path, "r");
    uint64_t total_size = 0;
    //指针定位到文件尾
    fseek(read_ref_file, 0, 2);
    //计算序列总长度
    total_size = ftell(read_ref_file);
    //指针定位到文件头
    fseek(read_ref_file, 0, 0);
    cout << "文件总长度是" << total_size << endl;

    char *p_ref_compress;
    p_ref_compress = (char*) malloc(sizeof(char) * total_size);
    fread(p_ref_compress, sizeof(char), total_size, read_ref_file);
    fclose(read_ref_file);

    //报存每一个$的位置
    map<uint64_t,uint64_t> my_count;

    uint64_t count = 0;
    uint64_t $_num = 0;
    my_count.insert(pair<uint64_t,uint64_t>(0,0));

    for(uint64_t i = 0; i < total_size; i++){
        if (p_ref_compress[i] == '\n' || p_ref_compress[i] == '\0') {
            continue;
        }
        if(p_ref_compress[i] == '$'){
            $_num++;
            my_count.insert(pair<uint64_t,uint64_t>(++count,i));
        }
    }
    cout << "$的数量是:" << $_num << endl;
    my_count.insert(pair<int,uint64_t>(++count,total_size-1));
    cout << "拼接前序列长度是：" << total_size << endl;
    return my_count;
}

