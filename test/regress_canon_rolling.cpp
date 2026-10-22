// Canonical rolling hashes for k > 32 must not depend on strand or on where a record starts.
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
static std::string revcomp(const std::string &s) {
    std::string ret(s.rbegin(), s.rend());
    for(auto &c: ret) c = c == 'A' ? 'T': c == 'C' ? 'G': c == 'G' ? 'C': 'A';
    return ret;
}
static std::set<uint64_t> canon_hashes(unsigned k, const std::string &s) {
    RollingHasher<uint64_t> rh(k, true);
    std::set<uint64_t> ret;
    rh.for_each_canon([&](uint64_t x) {ret.insert(x);}, s.data(), s.size());
    return ret;
}
TEST_CASE("A sequence and its reverse complement have the same canonical k-mers (k = 40)") {
    const std::string s = random_seq(1000, "ACGT", 1);
    REQUIRE(canon_hashes(40, s) == canon_hashes(40, revcomp(s)));
}
TEST_CASE("Dropping the first base keeps every later k-mer (k = 40)") {
    const std::string s = random_seq(1000, "ACGT", 2);
    const auto all = canon_hashes(40, s), rest = canon_hashes(40, s.substr(1));
    REQUIRE(std::includes(all.begin(), all.end(), rest.begin(), rest.end()));
}
