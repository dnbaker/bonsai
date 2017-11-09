#include "getopt.h"
#include "omp.h"
#include "lib/util.h"
#include <set>
#include <iostream>
#include <fstream>
#include <cassert>
#include <atomic>


#if 0

static std::atomic<uint64_t> n_inc = 0;
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
    std::cerr << "Incremented total of " << ++n_inc << "times. Expected " << nels << " times\n";
    return ret = map.n_occupied;
}

template<typename Map>
void reset(Map &map) {
    Map tmp;
    std::swap(tmp, map);
}

template<typename Map>
auto sorted_hash_els(const Map &map) {
    std::vector<int> ret;
    for(size_t i(0); i < map.size; ++i) {
        if(map.exist(i)) ret.push_back(map.vals[i]);
    }
    std::sort(std::begin(ret), std::end(ret));
    return ret;
}

int main(int argc, char *argv[]) {
    int nels(1 << 10), nthreads(1), c, csum(0), niter(100);
    emp::kh::khpp_t<int, int> map;
    static_assert(std::is_same<typename decltype(map)::val_type, int>::value, "Must have int  as a value.");
    std::fstream ofs("output.txt", std::ios_base::out);
    while((c = getopt(argc, argv, "i:n:p:h?")) >= 0) {
        switch(c) {
            case 'h': case '?': 
                std::fprintf(stderr, "-n:\t Number of elements. [1 << 10]\n-p:\t Number of processes. [1]\n");
                return EXIT_FAILURE;
            case 'n': nels = atoi(optarg); break;
            case 'p': nthreads = atoi(optarg); break;
            case 'i': niter = atoi(optarg); break;
        }
    }
    std::set<int> tels;
    for(int i(0); i < nels; ++i) tels.insert(i);
    std::vector<int> els;
    for(const auto i: tels) els.insert(els.begin(), i, i);
    std::cerr << "size of els " << els.size() << '\n';
    std::cerr << "size of tels " << tels.size() << '\n';
    std::random_shuffle(els.begin(), els.end());
    std::unordered_map<int, int> counts;
    int ret;
    omp_set_num_threads(nthreads);
    map = decltype(map)();
    size_t npadded;
    TIME_CODE(
    {
        for(int j(0); j < niter; ++j) {
            for(size_t i = 0; i < static_cast<unsigned>(nels); map.iput(els[i++], &ret));
            npadded = map.n_occupied;
        }
    }, "single build.");
    for(auto &el: els) {
        assert(map.iget(el) != map.n_buckets);
        std::cerr << "Count: " << map.vals[map.iget(el)]  << '\n';
        assert(map.vals[map.iget(el)] == 1);
    }
    reset(map);
    parallel_add_els(map, els, nthreads);
    npadded = map.n_occupied;
    auto sels(sorted_hash_els(map));
    std::cerr << "sorted in map: \n";
    for(auto el: sels) std::cerr << el << ", ";
    std::cerr << "\n";
    return EXIT_SUCCESS;
}
#else
int main() {
    std::cerr << "Hash map fails all tests. It does not work. Do not use. I will try to entirely rewrite it later.\n";
    return 1;
}
#endif
