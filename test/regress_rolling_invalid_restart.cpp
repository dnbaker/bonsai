#include "test/catch.hpp"
#include "encoder.h"
#include <random>
#include <set>
#include <string>
#include <vector>

using namespace bns;

// After an invalid character, the rolling hasher must resume with the k-mer
// that starts right after it: the hashes of a sequence with N characters are
// exactly the union of the hashes of its maximal N-free runs.

static std::string random_seq(size_t n, uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::string ret(n, 'A');
    for(auto &c: ret) c = "ACGT"[rng() % 4];
    return ret;
}
template<typename IntType>
static std::set<IntType> rolling_hashes(const std::string &s, unsigned k, bool canon) {
    RollingHasher<IntType> rh(k, canon);
    std::set<IntType> ret;
    rh.for_each_hash([&](IntType x) {ret.insert(x);}, s.data(), s.size());
    return ret;
}
template<typename IntType>
static std::set<IntType> union_of_runs(const std::string &s, unsigned k, bool canon) {
    std::set<IntType> ret;
    size_t start = 0;
    for(size_t i = 0; i <= s.size(); ++i) {
        if(i == s.size() || s[i] == 'N') {
            const auto h = rolling_hashes<IntType>(s.substr(start, i - start), k, canon);
            ret.insert(h.begin(), h.end());
            start = i + 1;
        }
    }
    return ret;
}
// Ns close together, far apart, and near the end of the sequence.
static std::string seq_with_ns(uint64_t seed) {
    std::string s = random_seq(600, seed);
    for(const size_t p: {150, 160, 300, 520}) s[p] = 'N';
    return s;
}

TEST_CASE("Rolling hashes after an N match the N-free runs (k = 40)") {
    const std::string s = seq_with_ns(1);
    for(const bool canon: {false, true}) {
        const auto expected = union_of_runs<uint64_t>(s, 40, canon);
        REQUIRE(expected.size() > 400);
        REQUIRE(rolling_hashes<uint64_t>(s, 40, canon) == expected);
    }
}
TEST_CASE("Rolling hashes after an N match the N-free runs (128-bit, k = 100)") {
    const std::string s = seq_with_ns(2);
    for(const bool canon: {false, true}) {
        const auto expected = union_of_runs<u128>(s, 100, canon);
        REQUIRE(expected.size() > 200);
        REQUIRE(rolling_hashes<u128>(s, 100, canon) == expected);
    }
}
