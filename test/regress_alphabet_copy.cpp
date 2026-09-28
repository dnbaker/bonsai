// Copies of an Encoder, and a RollingHasher built for protein, must parse with the protein table.
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
TEST_CASE("A copied PROTEIN20 encoder sees every 8-mer of a protein") {
    const std::string s = random_seq(500, "ACDEFGHIKLMNPQRSTVWY", 3);
    Encoder<score::Lex, uint64_t> orig(Spacer(8), false);
    orig.hashtype(PROTEIN20);
    Encoder<score::Lex, uint64_t> copied(orig);
    size_t n = 0;
    copied.for_each([&](uint64_t) {++n;}, s.data(), s.size());
    REQUIRE(n == s.size() - 8 + 1);
}
TEST_CASE("A PROTEIN20 RollingHasher sees every 20-mer of a protein") {
    const std::string s = random_seq(500, "ACDEFGHIKLMNPQRSTVWY", 4);
    RollingHasher<uint64_t> rh(20, false, PROTEIN20);
    size_t n = 0;
    rh.for_each_hash([&](uint64_t) {++n;}, s.data(), s.size());
    REQUIRE(n == s.size() - 20 + 1);
}
