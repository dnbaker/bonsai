#include "test/catch.hpp"
#include "encoder.h"
#include <algorithm>
#include <random>
#include <string>
#include <vector>

using namespace bns;

// Spaced (necessarily non-canonical) minimizers: a window covers the last
// w - c + 1 valid spaced k-mers (c the span), k-mers with an invalid
// character are skipped, a record with fewer valid k-mers than the window
// yields its best one, and poly-T (all ones when the k-mer fills the word)
// is emitted like any other k-mer.

static std::string random_seq(size_t n, uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::string ret(n, 'A');
    for(auto &c: ret) c = "ACGT"[rng() % 4];
    return ret;
}
static std::vector<unsigned> offsets(const std::vector<unsigned> &gaps) {
    std::vector<unsigned> ret{0};
    for(const auto g: gaps) ret.push_back(ret.back() + g + 1);
    return ret;
}
static std::string spacing_string(const std::vector<unsigned> &gaps) {
    std::string ret;
    for(const auto g: gaps) ret += (ret.empty() ? "": ",") + std::to_string(g);
    return ret;
}
// Forward 2-bit encodings of the valid spaced k-mers, in sequence order.
static std::vector<uint64_t> valid_kmers(const std::string &s, const std::vector<unsigned> &offs) {
    std::vector<uint64_t> ret;
    for(size_t i = 0; i + offs.back() < s.size(); ++i) {
        uint64_t f = 0;
        bool ok = true;
        for(const auto o: offs) {
            const auto c = std::string("ACGT").find(s[i + o]);
            if(c == std::string::npos) {ok = false; break;}
            f = (f << 2) | c;
        }
        if(ok) ret.push_back(f);
    }
    return ret;
}
static std::vector<uint64_t> expected_minimizers(const std::string &s, const std::vector<unsigned> &gaps, unsigned w) {
    const auto offs = offsets(gaps);
    const auto km = valid_kmers(s, offs);
    const size_t c = offs.back() + 1, m = std::max<size_t>(w, c) - c + 1;
    auto better = [](uint64_t x, uint64_t y) {return std::make_pair(lex_score(x), x) < std::make_pair(lex_score(y), y);};
    std::vector<uint64_t> ret;
    if(km.empty()) return ret;
    if(km.size() < m) return {*std::min_element(km.begin(), km.end(), better)};
    for(size_t j = m - 1; j < km.size(); ++j)
        ret.push_back(*std::min_element(km.begin() + (j + 1 - m), km.begin() + (j + 1), better));
    return ret;
}
// Runs the file API (as dashing2 does) on a one-record FASTA file.
static std::vector<uint64_t> minimizers(const std::string &s, const std::vector<unsigned> &gaps, unsigned w) {
    const unsigned k = gaps.size() + 1;
    const std::string path = "regress_spaced_window_invalid.tmp.fa";
    std::FILE *fp = std::fopen(path.data(), "w");
    std::fprintf(fp, ">r\n%s\n", s.data());
    std::fclose(fp);
    Encoder<score::Lex, uint64_t> enc(Spacer(k, w, spacing_string(gaps).data()), false);
    std::vector<uint64_t> ret;
    enc.for_each([&](uint64_t x) {ret.push_back(x);}, path.data());
    std::remove(path.data());
    return ret;
}

TEST_CASE("Spaced windows skip k-mers that contain N") {
    const std::vector<unsigned> gaps{0, 1, 0, 2, 0, 0, 1};
    std::string s = random_seq(400, 1);
    for(const size_t p: {60, 61, 62, 150, 151, 300}) s[p] = 'N';
    const auto got = minimizers(s, gaps, 20);
    const auto exp = expected_minimizers(s, gaps, 20);
    REQUIRE(std::find(got.begin(), got.end(), ~uint64_t(0)) == got.end());
    REQUIRE(got == exp);
}
TEST_CASE("A record shorter than a spaced window yields its best k-mer") {
    const std::vector<unsigned> gaps{0, 1, 0, 2, 0, 0, 1};
    const std::string s = random_seq(18, 2);  // 7 k-mers of span 12, window holds 9
    const auto exp = expected_minimizers(s, gaps, 20);
    REQUIRE(exp.size() == 1);
    REQUIRE(minimizers(s, gaps, 20) == exp);
}
TEST_CASE("Spaced poly-T k-mers that fill the word are emitted (k = 32)") {
    std::vector<unsigned> gaps(31, 0);
    gaps[10] = 1;
    const std::string s = random_seq(50, 3) + std::string(70, 'T') + random_seq(50, 4);
    for(const unsigned w: {0u, 40u}) {
        const auto got = minimizers(s, gaps, w);
        const auto exp = expected_minimizers(s, gaps, w);
        REQUIRE(std::count(exp.begin(), exp.end(), ~uint64_t(0)) > 0);
        REQUIRE(got == exp);
    }
}
