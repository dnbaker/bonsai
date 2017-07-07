#include "getopt.h"
#include "omp.h"
#include "lib/util.h"
#include "lib/khpp.h"
#include <iostream>

int main(int argc, char *argv[]) {
    int nels(1 << 10), nthreads(1), c;
    emp::kh::khpp_t<int, int> map;
    while((c = getopt(argc, argv, "n:p:h?")) >= 0) {
        switch(c) {
            case 'h': case '?': goto usage;
            case 'n': nels = atoi(optarg); break;
            case 'p': nthreads = atoi(optarg); break;
        }
    }
    #pragma omp parallel for num_threads(nthreads)
    for(int i = 0; i < nels; ++i) {
        int ret;
        std::cerr << "Index" << map.iput(i, &ret) << ".\n";
    }
    TIME_CODE(
    {
        auto it(map.begin());
        auto eit(map.end());
        while(it != eit) {
            std::cerr << "Key: " << it.first() << ". Value: " << it.second() << ".\n";
            ++it;
        }
    }
              , "khpp iteration.");
    TIME_CODE(
    {
        for(khiter_t ki = 0;ki != map.n_buckets; ++ki) if(map.exist(ki)) std::cerr << " Key: " << map.keys[ki] << ". Value: " << map.vals[ki] << ".\n";
    }
              , "kh iteration.");
    return EXIT_SUCCESS;
    usage:
    std::fprintf(stderr, "-n:\t Number of elements. [1 << 10]\n-p:\t Number of processes. [1]\n");
    return EXIT_FAILURE;
}
