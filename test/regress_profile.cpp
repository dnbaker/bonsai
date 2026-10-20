// Sample profiles (profile.h) and the counter state they rely on (uniq.h):
//  - UniqCounter::canonicalize: the state depends only on the set of hashes
//    (exact iff <= EXACT_MAX distinct), whatever the add order or duplicates;
//  - write/read round trip, exact and both sketch kinds, registers identical;
//  - counter_overlap: exact for exact counters, b-bit MinHash Jaccard within
//    0.03 of the truth at p = 12;
//  - abundance_distance: Bray-Curtis and tree EMD on a two-leaf tree.
#include "test/catch.hpp"

#include "profile.h"

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
