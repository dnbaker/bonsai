// Sample profiles (profile.h) and the counter state they rely on (uniq.h):
//  - UniqCounter::canonicalize: the state depends only on the set of hashes
//    (exact iff <= EXACT_MAX distinct), whatever the add order or duplicates;
//  - write/read round trip, exact and both sketch kinds, registers identical;
//  - counter_overlap: exact for exact counters, b-bit MinHash Jaccard within
//    0.03 of the truth at p = 12;
//  - abundance_distance: Bray-Curtis and tree EMD on a two-leaf tree.
#include "test/catch.hpp"

#include "profile_stats.h"

#include <cstdio>
#include <random>
#include <string>
#include <unistd.h>
#include <vector>

using namespace bns;

namespace {
std::vector<uint64_t> hashes(size_t n, uint64_t seed) {
    std::vector<uint64_t> v(n);
    for(auto &x: v) x = uniq_hash(seed++);
    return v;
}
std::string roundtrip_path() {return "/tmp/bns_regress_profile_" + std::to_string(getpid()) + ".gz";}
UniqCounter roundtrip(UniqCounter &c) {
    const std::string path = roundtrip_path();
    gzFile fp = gzopen(path.data(), "wb");
    c.write(fp);
    gzclose(fp);
    UniqCounter r;
    fp = gzopen(path.data(), "rb");
    r.read(fp);
    gzclose(fp);
    std::remove(path.data());
    return r;
}
}

TEST_CASE("UniqCounter: canonical state depends only on the hash set") {
    for(const UniqSketch kind: {UniqSketch::HLL, UniqSketch::BBMH}) {
        // 700 distinct: in order (no compaction past EXACT_MAX before 1024 adds)
        // vs each added three times shuffled (compactions while < EXACT_MAX).
        const auto h = hashes(700, 7);
        UniqCounter a(10, kind), b(10, kind);
        for(const auto x: h) a.add(x);
        std::vector<uint64_t> rep;
        for(int i = 0; i < 3; ++i) rep.insert(rep.end(), h.begin(), h.end());
        std::shuffle(rep.begin(), rep.end(), std::mt19937_64(3));
        for(const auto x: rep) b.add(x);
        a.canonicalize(); b.canonicalize();
        REQUIRE(a.sketched());
        REQUIRE(b.sketched());
        REQUIRE(a.estimate() == b.estimate());
        // At most EXACT_MAX distinct stays exact.
        UniqCounter e(10, kind);
        for(const auto x: hashes(UniqCounter::EXACT_MAX, 99)) e.add(x), e.add(x);
        REQUIRE(e.estimate() == double(UniqCounter::EXACT_MAX));
        REQUIRE(!e.sketched());
    }
}

TEST_CASE("UniqCounter: write/read round trip") {
    for(const UniqSketch kind: {UniqSketch::HLL, UniqSketch::BBMH}) {
        for(const size_t n: {size_t(0), size_t(300), size_t(5000)}) {
            UniqCounter c(11, kind);
            for(const auto x: hashes(n, 1)) c.add(x);
            UniqCounter r = roundtrip(c);
            REQUIRE(r.sketched() == c.sketched());
            REQUIRE(r.kind() == kind);
            REQUIRE(r.p() == 11);
            REQUIRE(r.estimate() == c.estimate());
            REQUIRE(r.exact() == c.exact());
            if(kind == UniqSketch::HLL && c.hll()) REQUIRE(r.hll()->core() == c.hll()->core());
            if(kind == UniqSketch::BBMH && c.bbmh_sketch()) REQUIRE(r.bbmh_sketch()->core() == c.bbmh_sketch()->core());
            // a loaded counter merges like a live one
            UniqCounter more(11, kind), live(11, kind);
            for(const auto x: hashes(2000, 50000)) more.add(x);
            for(const auto x: hashes(n, 1)) live.add(x);
            UniqCounter more2(11, kind);
            for(const auto x: hashes(2000, 50000)) more2.add(x);
            r.merge(more); live.merge(more2);
            REQUIRE(r.estimate() == live.estimate());
        }
    }
}

TEST_CASE("counter_overlap: exact and b-bit MinHash") {
    // A = [0, 400), B = [200, 600): exact, shared 200, Jaccard 1/3
    UniqCounter a(12, UniqSketch::BBMH), b(12, UniqSketch::BBMH);
    for(uint64_t i = 0; i < 400; ++i) a.add(uniq_hash(i));
    for(uint64_t i = 200; i < 600; ++i) b.add(uniq_hash(i));
    auto o = counter_overlap(a, b);
    REQUIRE(o.shared == 200.);
    REQUIRE(std::abs(o.jaccard - 1. / 3.) < 1e-12);
    // A = [0, 60000), B = [30000, 90000): Jaccard 1/3, shared 30000
    UniqCounter c(12, UniqSketch::BBMH), d(12, UniqSketch::BBMH);
    for(uint64_t i = 0; i < 60000; ++i) c.add(uniq_hash(i));
    for(uint64_t i = 30000; i < 90000; ++i) d.add(uniq_hash(i));
    o = counter_overlap(c, d);
    REQUIRE(std::abs(o.jaccard - 1. / 3.) < 0.03);
    REQUIRE(std::abs(o.shared / 30000. - 1.) < 0.1);
    // exact vs sketch: [0, 400) is contained in [0, 60000)
    o = counter_overlap(a, c);
    REQUIRE(std::abs(o.shared / 400. - 1.) < 0.25);
}

TEST_CASE("abundance_distance: Bray-Curtis and tree EMD") {
    // root 1 -> 2 -> {3, 4}; a: 10 reads at 3; b: 5 at 3, 5 at 4
    std::unordered_map<tax_t, uint64_t> a{{3, 10}}, b{{3, 5}, {4, 5}};
    std::unordered_map<tax_t, uint64_t> ca{{1, 10}, {2, 10}, {3, 10}}, cb{{1, 10}, {2, 10}, {3, 5}, {4, 5}};
    const auto d = abundance_distance(a, b, ca, 10, cb, 10, 1);
    REQUIRE(std::abs(d.bray_curtis - 0.5) < 1e-12);
    REQUIRE(std::abs(d.jaccard - 0.5) < 1e-12);
    REQUIRE(d.shared == 1);
    REQUIRE(std::abs(d.tree_emd - 1.) < 1e-12);   // half the mass moves 3 -> 2 -> 4
}

TEST_CASE("profile statistics: Clopper-Pearson, G, e-values, e-BH, lambda") {
    // Clopper-Pearson pinned to scipy.stats.beta.ppf (checked 2026-10-20)
    auto ci = clopper_pearson(3, 10);
    REQUIRE(std::abs(ci.first - 6.6739e-02) < 1e-6);
    REQUIRE(std::abs(ci.second - 6.5245e-01) < 1e-5);
    ci = clopper_pearson(0, 10);
    REQUIRE(ci.first == 0.);
    REQUIRE(std::abs(ci.second - 3.0850e-01) < 1e-5);
    ci = clopper_pearson(57, 400);
    REQUIRE(std::abs(ci.first - 1.0975e-01) < 1e-5);
    REQUIRE(std::abs(ci.second - 1.8065e-01) < 1e-5);
    // G statistic and its chi-square p (scipy chi2.sf)
    REQUIRE(std::abs(binomial_g(57, 400, 0.3) - 54.3512) < 1e-3);
    REQUIRE(std::abs(chisq1_p(binomial_g(10, 10, 0.3)) - 9.244e-07) < 1e-9);
    // e-value with rho = 0: the uniform-share Bayes factor 1 / ((n + 1) Binom(k; n, p0))
    {
        const double k = 7, n = 20, p0 = 0.2;
        const double lb = std::lgamma(n + 1) - std::lgamma(k + 1) - std::lgamma(n - k + 1) + k * std::log(p0) + (n - k) * std::log1p(-p0);
        REQUIRE(std::abs(log_evalue_share(k, n, p0, 0.) - (-std::log(n + 1) - lb)) < 1e-9);
        // with overdispersion the e-value shrinks (the null explains more)
        REQUIRE(log_evalue_share(k, n, p0, 0.05) < log_evalue_share(k, n, p0, 0.));
        // interval null containing the observed share: no evidence beyond the band's best null
        REQUIRE(log_evalue_share_band(10, 40, 0.2, 0.3, 0.) <= log_evalue_share(10, 40, 0.25, 0.) + 1e-9);
    }
    // e-BH: 4 e-values, alpha 0.05 -> reject the k largest with e_(k) >= m / (alpha k)
    REQUIRE(ebh_threshold({200, 50, 1, 0.5}, 0.05) == 4. / (0.05 * 2));   // e_(2) = 50 >= 40: reject 2
    REQUIRE(std::isinf(ebh_threshold({10, 5}, 0.05)));
    // lambda: hits = P lambda, distinct = P (1 - e^-lambda) -> recovered
    for(const double lam: {0.2, 1.0, 5.0}) {
        const double P = 1e5, hits = P * lam, distinct = P * -std::expm1(-lam);
        REQUIRE(std::abs(clade_lambda(hits, distinct) - lam) < 1e-6 * std::max(1., lam));
    }
    REQUIRE(std::isnan(clade_lambda(100, 100)));   // no duplicates: unidentifiable
}
