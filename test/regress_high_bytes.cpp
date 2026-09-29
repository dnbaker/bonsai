#include "test/catch.hpp"
#include "encoder.h"
#include <random>
#include <set>
#include <string>

using namespace bns;

// A byte in 0x80-0xff is not a residue, so it must end k-mers exactly as an
// invalid letter (X, which is neither a base nor an amino acid) does.
static std::string random_seq(size_t n, const std::string &alphabet, uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::string ret(n, 'A');
    for(auto &c: ret) c = alphabet[rng() % alphabet.size()];
    return ret;
}
template<typename KmerT>
static std::set<KmerT> encoder_kmers(const std::string &s, unsigned k, InputType rht, bool canon) {
    Encoder<score::Lex, KmerT> enc(Spacer(k), canon);
    enc.hashtype(rht);
    std::set<KmerT> ret;
    enc.for_each([&](KmerT x) {ret.insert(x);}, s.data(), s.size());
    return ret;
}
static std::set<uint64_t> rolling_hashes(const std::string &s, unsigned k, bool canon) {
    RollingHasher<uint64_t> rh(k, canon);
    std::set<uint64_t> ret;
    rh.for_each_hash([&](uint64_t x) {ret.insert(x);}, s.data(), s.size());
    return ret;
}
// Returns the high bytes that give different k-mers than an X in the same place.
template<typename F>
static std::string bytes_read_as_residues(std::string s, const F &kmers) {
    const size_t mid = s.size() / 2;
    s[mid] = 'X';
    const auto expected = kmers(s);
    std::string bad;
    for(unsigned b = 0x80; b < 0x100; ++b) {
        s[mid] = char(b);
        if(kmers(s) != expected) bad += std::to_string(b) + " ";
    }
    return bad;
}

TEST_CASE("High bytes are not read as bases by the DNA encoder (k = 21, and k = 40 in 128 bits)") {
    const std::string s = random_seq(301, "ACGT", 1);
    for(const bool canon: {false, true}) {
        REQUIRE(bytes_read_as_residues(s, [&](const std::string &x) {return encoder_kmers<uint64_t>(x, 21, DNA, canon);}) == "");
        REQUIRE(bytes_read_as_residues(s, [&](const std::string &x) {return encoder_kmers<u128>(x, 40, DNA, canon);}) == "");
    }
}
TEST_CASE("High bytes are not read as bases by the rolling hasher (k = 40)") {
    const std::string s = random_seq(301, "ACGT", 2);
    for(const bool canon: {false, true})
        REQUIRE(bytes_read_as_residues(s, [&](const std::string &x) {return rolling_hashes(x, 40, canon);}) == "");
}
TEST_CASE("High bytes are not read as residues by the protein encoder (k = 6)") {
    const std::string s = random_seq(301, "ACDEFGHIKLMNPQRSTVWY", 3);
    REQUIRE(bytes_read_as_residues(s, [&](const std::string &x) {return encoder_kmers<uint64_t>(x, 6, PROTEIN20, false);}) == "");
}
