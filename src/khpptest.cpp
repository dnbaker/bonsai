#include "getopt.h"
#include "omp.h"
#include "lib/util.h"
#include "lib/khpp.h"
#include <set>
#include <iostream>
#include <fstream>
#include <cassert>

int main(int argc, char *argv[]) {
    int nels(1 << 10), nthreads(1), c, csum(0);
    emp::kh::khpp_t<int, int> map;
    std::fstream ofs("output.txt", std::ios_base::out);
    while((c = getopt(argc, argv, "n:p:h?")) >= 0) {
        switch(c) {
            case 'h': case '?': 
                std::fprintf(stderr, "-n:\t Number of elements. [1 << 10]\n-p:\t Number of processes. [1]\n");
                return EXIT_FAILURE;
            case 'n': nels = atoi(optarg); break;
            case 'p': nthreads = atoi(optarg); break;
        }
    }
    std::set<int> els;
    while(els.size() < static_cast<unsigned>(nels)) els.insert(std::rand());
    int ret;
    omp_set_num_threads(nthreads);
    TIME_CODE(
    {
        auto it(els.begin());
        for(int i = 0; i < nels; ++i) {
            map.iput(*it++, &ret);
        }
    }, "parallel build.");
    map = decltype(map)();
    TIME_CODE(
    {
        auto it(els.begin());
        for(int i = 0; i < nels; ++i) {
            map.iput(*it++, &ret);
        }
    }, "single build.");
    TIME_CODE(
    {
        auto it(map.begin());
        auto eit(map.end());
        while(it != eit) {
            csum += it.second();
            ++it;
        }
    }
              , "khpp iteration.");
    TIME_CODE(
    {
        for(khiter_t ki = 0;ki != map.n_buckets; ++ki) if(map.exist(ki)) csum += map.vals[ki];
    }
              , "kh iteration.");
    std::set<int> cmp;
    for(auto it(map.begin()); it != map.end(); ++it) {
        cmp.insert(it.first());
    }
    assert(cmp == els);
    for(auto it(map.begin()), eit(map.end()); it != eit; ++it) ofs << "Outputting stuff. Key: " << it.first() << ". Value: " << it.second() << "\n";
    std::fprintf(stderr, "Sum of all the stuffs: %i\n", csum);
    return EXIT_SUCCESS;
}
