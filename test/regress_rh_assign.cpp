// Assigning a RollingHasher must adopt the assigned k.
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
TEST_CASE("RollingHasher assigned from k = 45 hashes like one constructed with k = 45") {
    const std::string s = random_seq(1000, "ACGT", 5);
    RollingHasher<uint64_t> assigned(21, true), fresh(45, true);
    assigned = fresh;
    std::set<uint64_t> a, b;
    assigned.for_each_canon([&](uint64_t x) {a.insert(x);}, s.data(), s.size());
    fresh.for_each_canon([&](uint64_t x) {b.insert(x);}, s.data(), s.size());
    REQUIRE(a == b);
}
