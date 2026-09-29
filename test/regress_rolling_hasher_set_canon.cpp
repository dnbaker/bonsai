// RollingHasherSet::for_each_canon must emit, for each of its k values, the same
// canonical hashes as a RollingHasher with that k and the same seeds, and so must
// not depend on strand.
#include "test/catch.hpp"
#include "encoder.h"
#include <map>
#include <random>
#include <string>
#include <vector>
using namespace bns;
static std::string random_dna(size_t n, uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::string ret(n, 'A');
    for(auto &c: ret) c = "ACGT"[rng() % 4];
    return ret;
}
static std::string revcomp(const std::string &s) {
    std::string ret(s.rbegin(), s.rend());
    for(auto &c: ret) c = c == 'A' ? 'T': c == 'C' ? 'G': c == 'G' ? 'C': 'A';
    return ret;
}
static const std::vector<unsigned> kvals{21, 33, 40};
static std::map<size_t, std::vector<uint64_t>> set_hashes(const std::string &s) {
    RollingHasherSet<uint64_t> rhs(kvals, true, DNA, 1337u);
    std::map<size_t, std::vector<uint64_t>> ret;
    rhs.for_each_canon([&](uint64_t x, size_t hi) {ret[hi].push_back(x);}, s.data(), s.size());
    return ret;
}
TEST_CASE("RollingHasherSet canonical hashes match RollingHasher for each k") {
    const std::string s = random_dna(1000, 7);
    const auto got = set_hashes(s);
    std::mt19937_64 mt(1337u);
    for(size_t hi = 0; hi < kvals.size(); ++hi) {
        const uint64_t seed1 = mt(), seed2 = mt();
        RollingHasher<uint64_t> rh(kvals[hi], true, DNA, -1, seed1, seed2);
        std::vector<uint64_t> expected;
        rh.for_each_canon([&](uint64_t x) {expected.push_back(x);}, s.data(), s.size());
        REQUIRE(expected.size() == s.size() - kvals[hi] + 1);
        REQUIRE(got.count(hi));
        REQUIRE(got.at(hi) == expected);
    }
}
TEST_CASE("RollingHasherSet canonical hashes do not depend on strand") {
    const std::string s = random_dna(1000, 8);
    const auto fwd = set_hashes(s), rev = set_hashes(revcomp(s));
    for(size_t hi = 0; hi < kvals.size(); ++hi) {
        REQUIRE(fwd.count(hi));
        REQUIRE(rev.count(hi));
        std::vector<uint64_t> r(rev.at(hi).rbegin(), rev.at(hi).rend());
        REQUIRE(fwd.at(hi) == r);
    }
}
