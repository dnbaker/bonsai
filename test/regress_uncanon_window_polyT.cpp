#include "test/catch.hpp"
#include "encoder.h"
#include <algorithm>
#include <random>
#include <string>
#include <vector>

using namespace bns;

// Non-canonical windowed minimizers when k fills the k-mer word (k = 32 in 64
// bits, 64 in 128 bits): poly-T then encodes as all ones, the value QueueMap
// also returns for a window that is not yet full. A full window whose
// minimizer is poly-T must still emit it.

static std::string random_seq(size_t n, uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::string ret(n, 'A');
    for(auto &c: ret) c = "ACGT"[rng() % 4];
    return ret;
}
// Forward 2-bit encodings of all k-mers, in sequence order (ACGT only).
template<typename KmerT>
static std::vector<KmerT> fwd_kmers(const std::string &s, unsigned k) {
    std::vector<KmerT> ret;
    for(size_t i = 0; i + k <= s.size(); ++i) {
        KmerT f = 0;
        for(size_t j = 0; j < k; ++j) f = (f << 2) | KmerT(std::string("ACGT").find(s[i + j]));
        ret.push_back(f);
    }
    return ret;
}
// One minimizer (lowest lex_score) per full window of w - k + 1 k-mers.
template<typename KmerT>
static std::vector<KmerT> expected_minimizers(const std::string &s, unsigned k, unsigned w) {
    const auto km = fwd_kmers<KmerT>(s, k);
    const size_t m = w - k + 1;
    std::vector<KmerT> ret;
    for(size_t j = m - 1; j < km.size(); ++j)
        ret.push_back(*std::min_element(km.begin() + (j + 1 - m), km.begin() + (j + 1),
            [](KmerT x, KmerT y) {return lex_score(x) < lex_score(y);}));
    return ret;
}
template<typename ScoreT, typename KmerT>
static std::vector<KmerT> minimizers(const std::string &s, unsigned k, unsigned w) {
    Encoder<ScoreT, KmerT> enc(Spacer(k, w), false);
    std::vector<KmerT> ret;
    enc.for_each([&](KmerT x) {ret.push_back(x);}, s.data(), s.size());
    return ret;
}

TEST_CASE("Non-canonical windows won by poly-T emit it (64-bit, k = 32)") {
    const std::string s = random_seq(50, 1) + std::string(60, 'T') + random_seq(50, 2);
    const auto got = minimizers<score::Lex, uint64_t>(s, 32, 40);
    const auto exp = expected_minimizers<uint64_t>(s, 32, 40);
    REQUIRE(std::count(exp.begin(), exp.end(), ~uint64_t(0)) > 0);
    REQUIRE(got == exp);
}
TEST_CASE("Non-canonical windows won by poly-T emit it (128-bit, k = 64)") {
    const std::string s = random_seq(50, 3) + std::string(90, 'T') + random_seq(50, 4);
    const auto got = minimizers<score::Lex, u128>(s, 64, 72);
    const auto exp = expected_minimizers<u128>(s, 64, 72);
    REQUIRE(std::count(exp.begin(), exp.end(), ~u128(0)) > 0);
    REQUIRE(got == exp);
}
TEST_CASE("Entropy-scored non-canonical windows emit one item per full window, poly-T included") {
    // Poly-T has zero entropy and the worst score, so it wins only the windows
    // that hold nothing else; there are 60 - 32 + 1 - 8 = 21 of them.
    const std::string s = random_seq(50, 5) + std::string(60, 'T') + random_seq(50, 6);
    const auto got = minimizers<score::Entropy, uint64_t>(s, 32, 40);
    REQUIRE(got.size() == s.size() - 40 + 1);
    REQUIRE(std::count(got.begin(), got.end(), ~uint64_t(0)) == 21);
}
