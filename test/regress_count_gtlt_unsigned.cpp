#include "test/catch.hpp"
#include "sketch/count_eq.h"
#include <random>

// count_gtlt on bytes and shorts compares the elements as unsigned values and
// must equal a plain loop. Values of 128 (32768 for shorts) and more are
// included. Build with the vector flags under test (for example -mavx512f
// without -mavx512bw).

template<typename T>
static size_t count_bad(std::mt19937_64 &rng) {
    alignas(64) static T a[4096 + 64], b[4096 + 64];
    size_t bad = 0;
    for(size_t n: {31, 64, 128, 200, 1000, 4096}) {
        for(size_t off: {0, 1}) {
            for(int rep = 0; rep < 20; ++rep) {
                for(size_t i = 0; i < n + off; ++i) {a[i] = T(rng()); b[i] = (rng() & 1) ? a[i]: T(rng());}
                uint64_t gt = 0, lt = 0;
                for(size_t i = 0; i < n; ++i) {gt += a[off + i] > b[off + i]; lt += a[off + i] < b[off + i];}
                const auto got = sketch::eq::count_gtlt(a + off, b + off, n);
                bad += got.first != gt || got.second != lt;
            }
        }
    }
    return bad;
}

TEST_CASE("count_gtlt on bytes and shorts matches a scalar count") {
    std::mt19937_64 rng(17);
    REQUIRE(count_bad<uint8_t>(rng) == 0);
    REQUIRE(count_bad<uint16_t>(rng) == 0);
}
