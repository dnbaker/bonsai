#include "test/catch.hpp"
#include "encoder.h"
#include <random>
#include <string>
#include <tuple>
#include <vector>

using namespace bns;

// Windowed canonical rolling hashing (k above the encoder's capacity) must
// select, for each window of w - k + 1 positions, the best of the per-position
// canonical hashes that the unwindowed hasher emits.

static std::string random_seq(size_t n, uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::string ret(n, 'A');
    for(auto &c: ret) c = "ACGT"[rng() % 4];
    return ret;
}
template<typename IntType>
static std::vector<IntType> hashes(unsigned k, long long w, const std::string &s) {
    RollingHasher<IntType> rh(k, true, DNA, w);
    std::vector<IntType> ret;
    rh.for_each_canon([&](IntType x) {ret.push_back(x);}, s.data(), s.size());
    return ret;
}
template<typename IntType>
static void check(unsigned k, unsigned w, uint64_t seed) {
    const std::string s = random_seq(3000, seed);
    const auto all = hashes<IntType>(k, -1, s);
    REQUIRE(all.size() == s.size() - k + 1);
    const size_t m = w - k + 1;
    std::vector<IntType> expected;
    for(size_t j = 0; j + m <= all.size(); ++j) {
        IntType best = all[j];
        for(size_t i = j + 1; i < j + m; ++i) {
            // Same order as the window: by score (lex_score, stored as uint64_t), then value.
            if(std::make_tuple(uint64_t(lex_score(all[i])), all[i]) < std::make_tuple(uint64_t(lex_score(best)), best))
                best = all[i];
        }
        expected.push_back(best);
    }
    REQUIRE(hashes<IntType>(k, w, s) == expected);
}

TEST_CASE("Windowed canonical rolling hashes: one canonical minimizer per window (64-bit, k = 40, w = 60)") {
    check<uint64_t>(40, 60, 1);
}
TEST_CASE("Windowed canonical rolling hashes: one canonical minimizer per window (128-bit, k = 100, w = 120)") {
    check<u128>(100, 120, 2);
}
