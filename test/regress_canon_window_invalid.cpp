#include "test/catch.hpp"
#include "encoder.h"
#include <algorithm>
#include <random>
#include <string>
#include <vector>

using namespace bns;

// Canonical windowed minimizers of a sequence containing N: k-mers that
// contain the N must not enter the window (they used to enter as the all-A
// k-mer, the canonical form of the overflow marker).

static std::string random_seq(size_t n, uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::string ret(n, 'A');
    for(auto &c: ret) c = "ACGT"[rng() % 4];
    return ret;
}
// Canonical 2-bit encodings of the k-mers without N, in sequence order.
template<typename KmerT>
static std::vector<KmerT> valid_canon_kmers(const std::string &s, unsigned k) {
    std::vector<KmerT> ret;
    for(size_t i = 0; i + k <= s.size(); ++i) {
        const std::string w = s.substr(i, k);
        if(w.find('N') != std::string::npos) continue;
        KmerT f = 0, r = 0;
        for(size_t j = 0; j < k; ++j) {
            const int c = std::string("ACGT").find(w[j]);
            f = (f << 2) | KmerT(c);
            r |= KmerT(3 - c) << (2 * j);
        }
        ret.push_back(std::min(f, r));
    }
    return ret;
}
// One minimizer (lowest lex_score) per full window of w - k + 1 valid k-mers.
template<typename KmerT>
static std::vector<KmerT> expected_minimizers(const std::string &s, unsigned k, unsigned w) {
    const auto km = valid_canon_kmers<KmerT>(s, k);
    const size_t m = w - k + 1;
    std::vector<KmerT> ret;
    for(size_t j = m - 1; j < km.size(); ++j)
        ret.push_back(*std::min_element(km.begin() + (j + 1 - m), km.begin() + (j + 1),
            [](KmerT x, KmerT y) {return lex_score(x) < lex_score(y);}));
    return ret;
}
template<typename KmerT>
static std::vector<KmerT> minimizers(const std::string &s, unsigned k, unsigned w) {
    Encoder<score::Lex, KmerT> enc(Spacer(k, w), true);
    std::vector<KmerT> ret;
    enc.for_each([&](KmerT x) {ret.push_back(x);}, s.data(), s.size());
    return ret;
}

TEST_CASE("Canonical windowed minimizers skip k-mers containing N (64-bit, k = 21 and 32)") {
    std::string s = random_seq(400, 1);
    s[200] = 'N';
    for(const unsigned k: {21u, 32u}) {
        const auto got = minimizers<uint64_t>(s, k, k + 5);
        REQUIRE(std::find(got.begin(), got.end(), uint64_t(0)) == got.end());
        REQUIRE(got == expected_minimizers<uint64_t>(s, k, k + 5));
    }
}
TEST_CASE("Canonical windowed minimizers skip k-mers containing N (128-bit, k = 40)") {
    std::string s = random_seq(400, 2);
    s[200] = 'N';
    const auto got = minimizers<u128>(s, 40, 45);
    REQUIRE(std::find(got.begin(), got.end(), u128(0)) == got.end());
    // Every minimizer is one of the sequence's canonical k-mers.
    Encoder<score::Lex, u128> enc(Spacer(40), true);
    std::vector<u128> all;
    enc.for_each([&](u128 x) {all.push_back(x);}, s.data(), s.size());
    for(const auto x: got) REQUIRE(std::find(all.begin(), all.end(), x) != all.end());
}
TEST_CASE("Poly-T still yields poly-A when the k-mer fills the word (k = 32)") {
    const std::string s = random_seq(20, 3) + std::string(40, 'T') + random_seq(20, 4);
    const auto got = minimizers<uint64_t>(s, 32, 34);
    REQUIRE(got == expected_minimizers<uint64_t>(s, 32, 34));
    REQUIRE(std::find(got.begin(), got.end(), uint64_t(0)) != got.end());
}
