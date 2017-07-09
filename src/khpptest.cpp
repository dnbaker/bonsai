#include "getopt.h"
#include "omp.h"
#include "lib/util.h"
#include "lib/khpp.h"
#include <set>
#include <iostream>
#include <fstream>
#include <cassert>

void parallel_add_els(emp::kh::khpp_t<int, int> &map, const std::vector<int> &els, int nthreads) {
    int ret;
    const size_t nels(els.size());
    #pragma omp parallel for num_threads(nthreads)
    for(size_t i = 0; i < nels; ++i) {
        map.iput(els[i], &ret);
    }
}

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
    std::set<int> tels;
    while(tels.size() < static_cast<unsigned>(nels)) tels.insert(std::rand());
    std::vector<int> els(tels.begin(), tels.end());
    int ret;
    omp_set_num_threads(nthreads);
    map = decltype(map)();
    TIME_CODE(
    {
        for(size_t i = 0; i < static_cast<unsigned>(nels); map.iput(els[i++], &ret));
    }, "single build.");
    map = decltype(map)();
    assert(map.size == 0);
    TIME_CODE(
    {
        parallel_add_els(map, els, nthreads);
    }, "parallel build.");

    TIME_CODE(
    {
        auto it(map.begin());
        auto eit(map.end());
        while(it != eit) {
            csum += it++.second();
        }
    }
              , "khpp iteration.");
    TIME_CODE(
    {
        for(khiter_t ki = 0;ki != map.n_buckets; ++ki) if(map.exist(ki)) csum += map.vals[ki];
    }
              , "kh iteration.");
    TIME_CODE(
    {
        for(const auto el: map) //std::cerr << "el.first: " << el.first << "el.second: " << el.second << ".\n";
            ret += el.first;
    }
              , "range for value");
    TIME_CODE(
    {
        for(const auto &el: map)
            // std::cerr << "el.first: " << el.first << "el.second: " << el.second << ".\n";
            ret += el.first;
    }
              , "range for reference");
    map = decltype(map)();
    els.clear();
    for(size_t i(0); i < nels; els.push_back(i++));
    for(size_t i = 0; i < nels; ++i) map.iput(i, &ret);
    std::vector<int> cmp;
    for(auto it(map.begin()); it != map.end(); ++it) {
        cmp.push_back(it.first());
    }
    if(cmp != els) {
        std::cerr << "cmp: \n";
        for(const auto i: cmp) std::cerr << i << ", ";
        std::cerr << '\n';
        std::cerr << "els: \n";
        for(const auto i: els) std::cerr << i << ", ";
        std::cerr << '\n';
    }
    for(auto it(map.begin()), eit(map.end()); it != eit; ++it) ofs << "Outputting stuff. Key: " << it.first() << ". Value: " << it.second() << "\n";
    map.clear();
    std::fprintf(stderr, "Sum of all the stuffs: %i\n", csum);
    return EXIT_SUCCESS;
}
