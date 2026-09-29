#include "test/catch.hpp"
#include "encoder.h"
#include <cmath>
#include <cstdio>
#include <random>
#include <string>
#include <unistd.h>
#include <vector>

using namespace bns;

// Entropy-weighted minimizers (score::Entropy) must prefer high-entropy
// k-mers: each window's minimizer minimizes hash / (entropy + 1e-4), where the
// hash is lex_score and the entropy is that of the k-mer's characters.

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
static uint64_t enc(const std::string &s) {
    uint64_t v = 0;
    for(char c: s) v = (v << 2) | std::string("ACGT").find(c);
    return v;
}
static double ref_score(const std::string &kmer) {
    double cnt[4] = {0, 0, 0, 0}, h = 0;
    for(char c: kmer) ++cnt[std::string("ACGT").find(c)];
    for(double x: cnt) if(x) h -= x / kmer.size() * std::log(x / kmer.size());
    return double(lex_score(enc(kmer)) >> 16) / (h + 1e-4);
}
// Checks that got[j] is a best-scoring k-mer of window j (up to rounding).
static void check_windows(const std::string &s, unsigned k, unsigned w, bool canon, const std::vector<uint64_t> &got) {
    const size_t m = w - k + 1, nk = s.size() - k + 1;
    REQUIRE(got.size() == nk - m + 1);
    for(size_t j = 0; j + m <= nk; ++j) {
        double best = INFINITY, gotscore = INFINITY;
        for(size_t i = j; i < j + m; ++i) {
            std::string km = s.substr(i, k);
            if(canon) km = std::min(km, revcomp(km));
            const double sc = ref_score(km);
            best = std::min(best, sc);
            if(enc(km) == got[j]) gotscore = sc;
        }
        REQUIRE(gotscore <= best * (1 + 1e-12) + 1);
    }
}
static std::string test_seq() {
    // Random sequence interleaved with low-entropy blocks.
    std::string s;
    for(int i = 0; i < 20; ++i)
        s += random_seq(40, i) + std::string(10 + i, "ACGT"[i % 4]) + std::string(i % 7 + 3, 'C') + "AT";
    return s;
}

TEST_CASE("String API, no canonicalization: minimizers minimize hash / entropy") {
    const std::string s = test_seq();
    Encoder<score::Entropy, uint64_t> e(Spacer(21, 42), false);
    std::vector<uint64_t> got;
    e.for_each([&](uint64_t x) {got.push_back(x);}, s.data(), s.size());
    check_windows(s, 21, 42, false, got);
}
TEST_CASE("File API: minimizers minimize hash / entropy, with and without canonicalization") {
    const std::string s = test_seq();
    char path[] = "/tmp/regress_entmin_XXXXXX";
    const int fd = mkstemp(path);
    REQUIRE(fd >= 0);
    const std::string rec = ">r\n" + s + "\n";
    REQUIRE(write(fd, rec.data(), rec.size()) == ssize_t(rec.size()));
    close(fd);
    for(const bool canon: {false, true}) {
        Encoder<score::Entropy, uint64_t> e(Spacer(21, 42), canon);
        std::vector<uint64_t> got;
        e.for_each([&](uint64_t x) {got.push_back(x);}, path);
        check_windows(s, 21, 42, canon, got);
    }
    std::remove(path);
}
TEST_CASE("A low-entropy k-mer is not chosen over high-entropy ones") {
    const std::string s = "AAAAAAAAAAAAAAAAAAAAC" "GTCAGTGCATCGATGCTAGCATGCAGTCGATGCACTG";
    Encoder<score::Entropy, uint64_t> e(Spacer(21, s.size()), false);
    std::vector<uint64_t> got;
    e.for_each([&](uint64_t x) {got.push_back(x);}, s.data(), s.size());
    REQUIRE(got.size() == 1);
    REQUIRE(got[0] != enc("AAAAAAAAAAAAAAAAAAAAC"));
}
