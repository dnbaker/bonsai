#include "diskarray.h"
using namespace emp;

int main(int argc, char *argv[]) {
    size_t n(argc == 1 ? 10000: std::atoi(argv[1]));
    ba::DiskBitArray ba(n, n, "whoo.bin");
    std::fprintf(stderr, "Tryna stuff\n");
    std::fprintf(stderr, "Popcount: %zu\n", ba.popcount());
    {
        Timer e(ks::sprintf("Set 1s for %zu by %zu", n, n).data());
        for(size_t i(0); i < n; ++i) {
            ba.set1(i);
        }
    }
    std::fprintf(stderr, "Popcount: %zu\n", ba.popcount());
}
