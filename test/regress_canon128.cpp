// 128-bit k-mers (33 <= k <= 64) must be canonicalized over their full length.
#include "test/catch.hpp"
#include "encoder.h"
#include <algorithm>
#include <random>
#include <set>
#include <string>
#include <vector>
using namespace bns;
static std::string random_seq(size_t n, const std::string &alphabet, uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::string ret(n, 'A');
    for(auto &c: ret) c = alphabet[rng() % alphabet.size()];
    return ret;
}
TEST_CASE("reverse_complement of a 40-mer of all A is a 40-mer of all T") {
    const u128 all_a = 0, all_t = (u128(1) << 80) - 1;  // A = 0, T = 3 in 2-bit encoding
    REQUIRE(reverse_complement(all_a, 40) == all_t);
    REQUIRE(canonical_representation(all_t, 40) == all_a);
}
TEST_CASE("k-mers that differ only in their first base canonicalize differently (k = 40)") {
    const u128 low = 0x0123456789abcdefULL;          // shared last 32 bases
    const u128 a = low, c = low | (u128(1) << 78);   // first base A vs C
    REQUIRE(canonical_representation(a, 40) != canonical_representation(c, 40));
}
