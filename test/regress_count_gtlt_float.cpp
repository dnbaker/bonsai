#include "test/catch.hpp"
#include "sketch/count_eq.h"
#include <random>

// count_gtlt on floats must equal a plain loop for every length, including
// lengths that leave whole vector registers after the unrolled loop and a
// remainder of fewer than one register. Build with the vector flags under
// test (for example -mavx2 without AVX-512).

TEST_CASE("count_gtlt on floats matches a scalar count") {
    std::mt19937_64 rng(19);
    alignas(64) static float a[400], b[400];
    size_t bad = 0;
    for(size_t n = 1; n < 300; ++n) {
        for(size_t i = 0; i < n; ++i) {a[i] = float(rng() % 1000); b[i] = (rng() & 1) ? a[i]: float(rng() % 1000);}
        uint64_t gt = 0, lt = 0;
        for(size_t i = 0; i < n; ++i) {gt += a[i] > b[i]; lt += a[i] < b[i];}
        const auto got = sketch::eq::count_gtlt(a, b, n);
        bad += got.first != gt || got.second != lt;
    }
    INFO(bad << " of 299 lengths differ");
    REQUIRE(bad == 0);
}
