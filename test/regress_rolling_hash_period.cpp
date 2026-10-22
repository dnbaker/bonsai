#include "test/catch.hpp"
#include "encoder.h"
#include <random>
#include <set>
#include <string>
#include <vector>

using namespace bns;

// Rolling hashes for k above the exact-encoding limit must not collide
// structurally for k larger than the hash word size.

static std::string random_seq(size_t n, uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::string ret(n, 'A');
    for(auto &c: ret) c = "ACGT"[rng() % 4];
    return ret;
}
template<typename IntType>
static std::set<IntType> hashes(unsigned k, const std::string &s) {
    RollingHasher<IntType> rh(k, /*canon=*/false);
    std::set<IntType> ret;
    rh.for_each_uncanon([&](IntType x) {ret.insert(x);}, s.data(), s.size());
    return ret;
}
template<typename IntType>
static void check_period(unsigned period, unsigned k) {
    std::string a, b;
    const std::string pa = random_seq(period, 1), pb = random_seq(period, 2);
    for(int i = 0; i < 8; ++i) a += pa, b += pb;
    const auto ha = hashes<IntType>(k, a), hb = hashes<IntType>(k, b);
    // A period-p sequence has exactly p distinct k-mers, one per phase.
    REQUIRE(ha.size() == period);
    REQUIRE(hb.size() == period);
    for(const auto x: ha) REQUIRE(hb.count(x) == 0);
}

TEST_CASE("Unrelated period-64 repeats do not collide (64-bit, k = 128)") {
    check_period<uint64_t>(64, 128);
}
TEST_CASE("Unrelated period-128 repeats do not collide (128-bit, k = 256)") {
    check_period<u128>(128, 256);
}
TEST_CASE("Swapping positions 0 and 64 changes the hash (k = 100)") {
    std::string a = random_seq(100, 3);
    if(a[0] == a[64]) a[0] = a[64] == 'A' ? 'C': 'A';
    std::string b = a;
    std::swap(b[0], b[64]);
    REQUIRE(hashes<uint64_t>(100, a) != hashes<uint64_t>(100, b));
    REQUIRE(hashes<u128>(100, a) != hashes<u128>(100, b));
}
TEST_CASE("Rolling updates agree with hashing each k-mer from scratch") {
    const unsigned k = 70;
    const std::string s = random_seq(400, 4);
    PolyHash<u128> fwd(k), rev(k), fresh(k);
    fwd.seed(1337, 137); rev.seed(1337, 137); fresh.seed(1337, 137);
    auto rc = [](char c) {return uint8_t(3 - std::string("ACGT").find(c));};
    auto code = [](char c) {return uint8_t(std::string("ACGT").find(c));};
    for(unsigned i = 0; i < k; ++i) fwd.eat(code(s[i]));
    for(unsigned i = k; i-- > 0;) rev.eat(rc(s[i]));
    for(size_t i = k; i <= s.size(); ++i) {
        fresh.reset();
        for(size_t j = i - k; j < i; ++j) fresh.eat(code(s[j]));
        REQUIRE(fwd.hashvalue == fresh.hashvalue);
        fresh.reset();
        for(size_t j = i; j-- > i - k;) fresh.eat(rc(s[j]));
        REQUIRE(rev.hashvalue == fresh.hashvalue);
        if(i == s.size()) break;
        fwd.update(code(s[i - k]), code(s[i]));
        rev.reverse_update(rc(s[i]), rc(s[i - k]));
    }
}
