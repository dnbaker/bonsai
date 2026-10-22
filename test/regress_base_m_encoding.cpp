// 20-, 14- and 6-letter alphabets are base-m integers: appending must add, not OR, and must not overflow.
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
static std::vector<u128> kmers(const std::string &s, unsigned k, InputType rht) {
    std::vector<u128> ret;
    if(k <= rh2n(rht, 8)) {
        Encoder<score::Lex, uint64_t> enc(Spacer(k), false); enc.hashtype(rht);
        enc.for_each([&](uint64_t x) {ret.push_back(x);}, s.data(), s.size());
    } else {
        Encoder<score::Lex, u128> enc(Spacer(k), false); enc.hashtype(rht);
        enc.for_each([&](u128 x) {ret.push_back(x);}, s.data(), s.size());
    }
    return ret;
}
static std::vector<u128> reference(const std::string &s, unsigned k, InputType rht) {
    std::vector<u128> ret;
    for(size_t i = 0; i + k <= s.size(); ++i) {
        u128 v = 0;
        for(size_t j = i; j < i + k; ++j) v = v * mul(rht) + u128(rh2lp(rht)[uint8_t(s[j])]);
        ret.push_back(v);
    }
    return ret;
}
TEST_CASE("PROTEIN20: CA and CF encode differently (20 | 4 == 20)") {
    REQUIRE(kmers("CA", 2, PROTEIN20) != kmers("CF", 2, PROTEIN20));
}
TEST_CASE("PROTEIN20 14-mers (the 64-bit limit) match base-20 arithmetic") {
    const std::string s = random_seq(300, "ACDEFGHIKLMNPQRSTVWY", 6);
    REQUIRE(kmers(s, 14, PROTEIN20) == reference(s, 14, PROTEIN20));
}
TEST_CASE("PROTEIN20 29-mers (the 128-bit limit) match base-20 arithmetic") {
    const std::string s = random_seq(300, "ACDEFGHIKLMNPQRSTVWY", 7);
    REQUIRE(kmers(s, 29, PROTEIN20) == reference(s, 29, PROTEIN20));
}
TEST_CASE("rhmask(PROTEIN20, 23) is exactly 20^23") {
    u128 p = 1; for(int i = 0; i < 23; ++i) p *= 20;
    REQUIRE(rhmask<u128>(PROTEIN20, 23) == p);
}
TEST_CASE("Encoder::kmer() on PROTEIN20 matches base-20 arithmetic") {
    const std::string s = "ACDEF";
    Encoder<score::Lex, uint64_t> enc(Spacer(5), false); enc.hashtype(PROTEIN20);
    enc.assign(s.data(), s.size());
    REQUIRE(u128(enc.kmer(0)) == reference(s, 5, PROTEIN20)[0]);
}
