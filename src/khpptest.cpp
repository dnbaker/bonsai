#include "getopt.h"
#include "omp.h"
#include "lib/util.h"
#include "lib/khpp.h"
#include <set>
#include <iostream>
#include <fstream>
#include <cassert>
template<typename T> class TD;

void parallel_add_els(emp::kh::khpp_t<int, int> &map, const std::vector<int> &els, int nthreads) {
    using map_type = emp::kh::khpp_t<int, int>;
    int ret;
    const size_t nels(els.size());
    nthreads = 1;
    for(size_t i = 0; i < nels; ++i) {
        auto lambda = [](const typename map_type::key_type &key, typename map_type::val_type &el){
            std::cerr << "Trying to execute.\n";
            __sync_fetch_and_add(&el, 1);
        };
        std::cerr << "made lambda\n";
        map.upsert(els[i], lambda, 0);
        assert(map.iget(els[i]) != map.n_buckets);
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
    std::cerr << std::boolalpha;
    omp_set_num_threads(nthreads);
    map = decltype(map)();
    TIME_CODE(
    {
        for(int i = 0; i < static_cast<unsigned>(nels); map.iput(els[i++], &ret))
            std::cerr << "Inserting element " << i << ".\n";
    }, "single build.");
    map = decltype(map)();
    for(const auto i: tels) {
        if(map.iget(i) == map.n_buckets) {
            for(khiter_t ki(0); ki < map.n_buckets; ++ki) {
                if(map.keys[ki] == i) {
                    std::cerr << "Found key with is del ? " << !!__ac_isdel(map.flags, ki) << " is empty " << !!__ac_isempty(map.flags, ki) << ".\n";
                    break;
                }
                std::cerr << "COuld not find key even with wrong flags.\n";
            }
        }
    }
    for(auto i: map) {
        std::cerr << "Find element " << i.first << " in map.\n";
        assert(tels.find(i.first) != tels.end());
    }
    std::cerr << "Single build worked.\n";
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
    for(auto &pair: counts) {
        auto m(map.iget(pair.first));
        std::cerr << "Key: " << pair.first << ". Count in 1 thread: " << pair.second << ". Count in parallel: " << (m == map.n_buckets ? std::string("missing"): std::to_string(m)) << ".\n";
    }
    return EXIT_SUCCESS;
}
