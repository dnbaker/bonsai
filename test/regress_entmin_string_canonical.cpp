#include "test/catch.hpp"
#include "encoder.h"
#include <random>
#include <set>
#include <string>

using namespace bns;

// Canonical entropy-weighted minimizers from Encoder::for_each(func, str, len)
// must not depend on the strand of the input.

static std::string random_seq(size_t n, uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::string ret(n, 'A');
    for(auto &c: ret) c = "ACGT"[rng() % 4];
    return ret;
}
static std::string revcomp(const std::string &s) {
    std::string ret(s.rbegin(), s.rend());
    for(auto &c: ret) c = c == 'A' ? 'T': c == 'C' ? 'G': c == 'G' ? 'C': 'A';
    return ret;
}
template<typename KmerT>
static std::set<KmerT> entmin_kmers(unsigned k, unsigned w, const std::string &s) {
    Encoder<score::Entropy, KmerT> enc(Spacer(k, w), /*canonicalize=*/true);
    std::set<KmerT> ret;
    enc.for_each([&](KmerT x) {ret.insert(x);}, s.data(), s.size());
    return ret;
}

TEST_CASE("String API entropy minimizers are strand-independent (64-bit, k = 21, w = 27)") {
    const std::string s = random_seq(3000, 21);
    const auto fwd = entmin_kmers<uint64_t>(21, 27, s);
    REQUIRE(fwd.size() > 100);
    REQUIRE(fwd == entmin_kmers<uint64_t>(21, 27, revcomp(s)));
}
