#include "test/catch.hpp"
#include "encoder.h"
#include <random>
#include <set>
#include <string>

using namespace bns;

// Encoder::for_each(func, str, len) must emit spaced k-mers, matching a
// direct encoding of every placement of the seed.

static std::string random_seq(size_t n, uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::string ret(n, 'A');
    for(auto &c: ret) c = "ACGT"[rng() % 4];
    return ret;
}

TEST_CASE("String API emits spaced k-mers (k = 12, seed 0,1,0,2,0x7)") {
    const std::string s = random_seq(500, 13);
    const spvec_t gaps{0, 1, 0, 2, 0, 0, 0, 0, 0, 0, 0};
    // Reference: 2-bit codes (A=0, C=1, G=2, T=3) at the kept positions,
    // first position most significant.
    std::set<uint64_t> expected;
    const size_t span = comb_size(gaps);
    for(size_t i = 0; i + span <= s.size(); ++i) {
        size_t pos = i;
        uint64_t v = std::string("ACGT").find(s[pos]);
        for(auto g: gaps) {
            pos += g + 1;
            v = (v << 2) | std::string("ACGT").find(s[pos]);
        }
        expected.insert(v);
    }
    Encoder<> enc(Spacer(12, 12, gaps), /*canonicalize=*/false);
    std::set<uint64_t> got;
    enc.for_each([&](uint64_t x) {got.insert(x);}, s.data(), s.size());
    REQUIRE(expected.size() > 400);
    REQUIRE(got == expected);
}
