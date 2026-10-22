#include "test/catch.hpp"
#include "encoder.h"
#include <random>
#include <set>
#include <string>

using namespace bns;

// Alphabet aliases: U reads as T in DNA (RNA input), and O and U read as K and C
// in protein alphabets (pyrrolysine and selenocysteine).

static std::string random_seq(size_t n, const std::string &alphabet, uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::string ret(n, 'A');
    for(auto &c: ret) c = alphabet[rng() % alphabet.size()];
    return ret;
}
static std::string replace_all(std::string s, char from, char to) {
    for(auto &c: s) if(c == from) c = to;
    return s;
}
static std::set<uint64_t> encoder_kmers(const std::string &s, unsigned k, InputType rht, bool canon) {
    Encoder<score::Lex, uint64_t> enc(Spacer(k), canon);
    enc.hashtype(rht);
    std::set<uint64_t> ret;
    enc.for_each([&](uint64_t x) {ret.insert(x);}, s.data(), s.size());
    return ret;
}
static std::set<uint64_t> rolling_hashes(const std::string &s, unsigned k, InputType rht, bool canon) {
    RollingHasher<uint64_t> rh(k, canon, rht);
    std::set<uint64_t> ret;
    rh.for_each_hash([&](uint64_t x) {ret.insert(x);}, s.data(), s.size());
    return ret;
}

TEST_CASE("RNA U and u encode as T (k = 21)") {
    const std::string dna = random_seq(300, "ACGT", 1);
    const std::string rna = replace_all(dna, 'T', 'U'), lower = replace_all(dna, 'T', 'u');
    for(const bool canon: {false, true}) {
        const auto expected = encoder_kmers(dna, 21, DNA, canon);
        REQUIRE(expected.size() > 250);
        REQUIRE(encoder_kmers(rna, 21, DNA, canon) == expected);
        REQUIRE(encoder_kmers(lower, 21, DNA, canon) == expected);
    }
}
TEST_CASE("RNA U hashes as T in the rolling hasher (k = 40)") {
    const std::string dna = random_seq(300, "ACGT", 2);
    const std::string rna = replace_all(dna, 'T', 'U'), lower = replace_all(dna, 'T', 'u');
    for(const bool canon: {false, true}) {
        const auto expected = rolling_hashes(dna, 40, DNA, canon);
        REQUIRE(expected.size() > 250);
        REQUIRE(rolling_hashes(rna, 40, DNA, canon) == expected);
        REQUIRE(rolling_hashes(lower, 40, DNA, canon) == expected);
    }
}
TEST_CASE("Protein O encodes as K and U as C") {
    const std::string p = random_seq(300, "ACDEFGHIKLMNPQRSTVWY", 3);
    const std::string po = replace_all(p, 'K', 'O'), pu = replace_all(p, 'C', 'U');
    for(const InputType rht: {PROTEIN20, PROTEIN_14, PROTEIN_6}) {
        const auto expected = encoder_kmers(p, 5, rht, false);
        REQUIRE(expected.size() > 100);
        REQUIRE(encoder_kmers(po, 5, rht, false) == expected);
        REQUIRE(encoder_kmers(pu, 5, rht, false) == expected);
    }
}
