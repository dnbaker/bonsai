#include "test/catch.hpp"
#include "sketch/count_eq.h"
#include <random>
#include <vector>

// count_eq on bytes must equal a plain loop over the elements. The vector
// branches process 32 or 64 bytes at a time and unroll by four, so lengths of
// 128 bytes and more reach the unrolled loop; both aligned and unaligned
// arrays are used. Build with the vector flags under test (for example -mavx2).

TEST_CASE("count_eq on bytes matches a scalar count") {
    std::mt19937_64 rng(13);
    alignas(64) static uint8_t a[4096 + 64], b[4096 + 64];
    size_t bad = 0, trials = 0;
    for(size_t n: {31, 32, 127, 128, 200, 1000, 1024, 4096}) {
        for(size_t off: {0, 1}) {
            for(int rep = 0; rep < 20; ++rep) {
                for(size_t i = 0; i < n + off; ++i) {a[i] = rng() & 3; b[i] = (rng() & 1) ? a[i]: rng() & 3;}
                size_t want = 0;
                for(size_t i = 0; i < n; ++i) want += a[off + i] == b[off + i];
                bad += sketch::eq::count_eq(a + off, b + off, n) != want;
                ++trials;
            }
        }
    }
    INFO(bad << " of " << trials << " counts differ");
    REQUIRE(bad == 0);
}
