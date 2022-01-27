//
// Created by lin on 21-10-25.
//
#include "basic.h"
#include "compress.h"
#include "split.h"
#include "utils.h"
#include "FM-index.h"

int main(int argc, char **argv) {

    //开始时间
    struct timeval tvs, tve;
    gettimeofday(&tvs, NULL);
    cout << "start..." << endl;
    if (strcmp(argv[1], "-c") == 0) {
        //执行NREPA算法即非重复边路径的抽取
        compress_main(argc, argv);
    }
    if (strcmp(argv[1], "-s") == 0) {
        //执行NREPS算法即非重复边路径的拼接
        split_main(argc,argv);
    }
    if (strcmp(argv[1], "-u") == 0) {
        //编写一些使用到的工具
        utils_main(argc, argv);
    }
    if (strcmp(argv[1], "-f") == 0) {
        //编写一些使用到的工具
        fm_main(argc, argv);
    }
    if (strcmp(argv[1], "-d") == 0) {
        //编写一些使用到的工具
        fm_main(argc, argv);
    }

    //结束时间
    cout << "end..." << endl;
    gettimeofday(&tve, NULL);
    double span = tve.tv_sec - tvs.tv_sec
                  + (tve.tv_usec - tvs.tv_usec) / 1000000.0;
    cout << "end time is: " << span << endl;
}