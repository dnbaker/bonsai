#include "test/catch.hpp"
#include "encoder.h"
#include <algorithm>
#include <random>
#include <string>
#include <vector>

using namespace bns;

// A record with fewer k-mers than the window yields one minimizer, the
// best-scoring k-mer of the record, whether or not k-mers are canonicalized.

static std::string random_seq(size_t n, uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::string ret(n, 'A');
    for(auto &c: ret) c = "ACGT"[rng() % 4];
    return ret;
}
template<typename KmerT>
static std::vector<KmerT> minimizers(const std::string &s, unsigned k, unsigned w, bool canon) {
    Encoder<score::Lex, KmerT> enc(Spacer(k, w), canon);
    std::vector<KmerT> ret;
    enc.for_each([&](KmerT x) {ret.push_back(x);}, s.data(), s.size());
    return ret;
}
template<typename KmerT>
static KmerT best_kmer(const std::string &s, unsigned k, bool canon) {
    Encoder<score::Lex, KmerT> enc(Spacer(k), canon);
    std::vector<KmerT> all;
    enc.for_each([&](KmerT x) {all.push_back(x);}, s.data(), s.size());
    return *std::min_element(all.begin(), all.end(), [](KmerT x, KmerT y) {return lex_score(x) < lex_score(y);});
}

TEST_CASE("A 40 bp record with k = 21 and w = 100 yields one minimizer (64-bit)") {
    const std::string s = random_seq(40, 1);
    for(const bool canon: {false, true}) {
        const auto got = minimizers<uint64_t>(s, 21, 100, canon);
        REQUIRE(got.size() == 1);
        REQUIRE(got[0] == best_kmer<uint64_t>(s, 21, canon));
    }
}
TEST_CASE("A 50 bp record with k = 40 and w = 100 yields one minimizer (128-bit)") {
    const std::string s = random_seq(50, 2);
    for(const bool canon: {false, true}) {
        const auto got = minimizers<u128>(s, 40, 100, canon);
        REQUIRE(got.size() == 1);
        REQUIRE(got[0] == best_kmer<u128>(s, 40, canon));
    }
}
