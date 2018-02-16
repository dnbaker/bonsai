#include "bitmap.h"
using namespace emp;
using namespace popcnt;

template<typename T>
unsigned vpopcnt(const T &con) {
    return vec_popcnt(con);
}

template<typename T>
size_t aln(T thing) {
    size_t ret(sizeof(T) * 8 - 1);
    uint64_t val((uint64_t) thing);
    while(val & ((1ULL << ret) - 1)) val >>= 1, --ret;
    return ret;
}

#define PELS std::fprintf(stderr, "Popcount is now %u\n", vpopcnt(els));

int main() {
    //lazy::vector<u64> els{0,0,0};
    lazy::vector<u64> els;
    while(els.size() < 3) els.push_back(0);
    for(auto el: els) {
        std::fprintf(stderr, "el: %zu\n", el);
    }
    PELS
    //els[2] |= 3;
    PELS
    els[0] |= ((1 << 16) | (1 << 8));
    PELS
    unsigned i(0);
    for(auto &el: els) i += __builtin_popcountll(el);
    std::fprintf(stderr, "Actual summed popcount: %u\n", i);
    auto data = new int[4000];
    std::fprintf(stderr, "Alignment: %zu. Number of empty bits: %u. %" PRIx64 "\n", aln(data), __builtin_ctz((uint64_t)data), (uint64_t)data);
    for(size_t i(0); i < 4000; ++i) data[i] = i;
    delete[] data;
    return 0;
}
