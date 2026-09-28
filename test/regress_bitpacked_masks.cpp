// Bit-packed alphabets must keep every character of a k-mer.
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
static std::set<uint64_t> kmers(const std::string &s, unsigned k, InputType rht) {
    Encoder<score::Lex, uint64_t> enc(Spacer(k), false);
    enc.hashtype(rht);
    std::set<uint64_t> ret;
    enc.for_each([&](uint64_t x) {ret.insert(x);}, s.data(), s.size());
    return ret;
}
TEST_CASE("PROTEIN_3BIT: 21-mers differing only in the first residue encode differently") {
    const std::string tail = "DEFGHIKLMNPQRSTVWYAC";
    REQUIRE(kmers("A" + tail, 21, PROTEIN_3BIT) != kmers("C" + tail, 21, PROTEIN_3BIT));
}
TEST_CASE("PROTEIN_3BIT: a 64-bit word holds at most 21 residues") {
    REQUIRE(rh2n(PROTEIN_3BIT, 8) * 3 <= 64);
}
TEST_CASE("DNAC: CCCCCCCCCC and AAAAAAAAAA encode differently") {
    REQUIRE(kmers("CCCCCCCCCC", 10, DNAC) != kmers("AAAAAAAAAA", 10, DNAC));
}
