#include "getopt.h"
#include "omp.h"
#include "lib/util.h"
#include "lib/khpp.h"
#include <set>
#include <iostream>
#include <fstream>
#include <cassert>

auto parallel_add_els(emp::kh::khpp_t<int, int> &map, const std::vector<int> &els, int nthreads) {
    using map_type = emp::kh::khpp_t<int, int>;
    int ret;
    const size_t nels(els.size());
    #pragma omp parallel num_threads(nthreads)
    for(size_t i = 0; i < nels; ++i) {
        map.upsert(els[i], [](const typename map_type::key_type &key, typename map_type::val_type &el){
            std::cerr << "Trying to execute. Old key: " << key << ". Old val: " << el << '\n';
            __sync_fetch_and_add(&el, 1);
            std::cerr << "Tried to execute. New key: " << key << ". New val: " << el << '\n';
        }, 0);
    }
    return ret = map.n_occupied;
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
    for(int i(0); i < nels; ++i) tels.insert(i);
    std::vector<int> els;
    for(const auto i: tels) els.insert(els.begin(), std::rand() & 0xF, i);
    std::cerr << "size of els " << els.size() << '\n';
    std::cerr << "size of tels " << tels.size() << '\n';
    std::random_shuffle(els.begin(), els.end());
    std::unordered_map<int, int> counts;
    for(auto i: els) ++counts[i];
    for(auto &pair: counts) std::cerr << pair.first << ", " << pair.second << '\n';
    int ret;
    omp_set_num_threads(nthreads);
    map = decltype(map)();
    size_t npadded;
    TIME_CODE(
    {
        for(size_t i = 0; i < static_cast<unsigned>(nels); map.iput(els[i++], &ret));
        npadded = map.n_occupied;
    }, "single build.");
    std::cerr << "single num added." << npadded << "\n";
    map = decltype(map)();
    assert(map.size == 0);
    TIME_CODE(
    {
        npadded = parallel_add_els(map, els, nthreads);
    }, "parallel build.");
    std::cerr << "parallel num added." << npadded << "\n";
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
    for(int i(0); i < nels; els.push_back(i++));
    for(int i(0); i < nels; map.iput(i++, &ret));
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
    std::fprintf(stderr, "Sum of all the stuffs: %i\n", csum);
    for(auto &pair: counts) {
        auto m(map.iget(pair.first));
        std::cerr << "Key: " << pair.first << ". Count in 1 thread: " << pair.second << ". Count in parallel: " << (m == map.n_buckets ? std::string("missing"): std::to_string(m)) << ".\n";
    }
    for(khiter_t i(0); i < map.n_buckets; ++i) if(map.exist(i)) std::cerr << "Key: " << map.keys[i] << " Val: " << map.vals[i] << "\n.";
    map.clear();
    return EXIT_SUCCESS;
}
