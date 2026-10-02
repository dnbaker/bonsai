#include "test/catch.hpp"
#include "encoder.h"
#include <random>
#include <string>
#include <vector>

using namespace bns;

// A windowed RollingHasher must not carry its window from one record to the
// next. Callers such as dashing2's --parse-by-seq read the window after a
// record that emitted nothing (one shorter than k) to pick a fallback item.

static std::string random_seq(size_t n, uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::string ret(n, 'A');
    for(auto &c: ret) c = "ACGT"[rng() % 4];
    return ret;
}

TEST_CASE("Non-canonical windowed RollingHasher starts each record with an empty window") {
    const std::string longrec = random_seq(200, 1), shortrec = random_seq(10, 2);
    for(const bool canon: {false, true}) {
        RollingHasher<uint64_t> rh(40, canon, DNA, 50);
        size_t n = 0;
        rh.for_each_hash([&](uint64_t) {++n;}, longrec.data(), longrec.size());
        if(!canon) REQUIRE(n == 200 - 50 + 1);
        n = 0;
        rh.for_each_hash([&](uint64_t) {++n;}, shortrec.data(), shortrec.size());
        REQUIRE(n == 0);
        REQUIRE(rh.n_in_queue() == 0);
    }
}
TEST_CASE("128-bit and protein RollingHasher start each record with an empty window") {
    const std::string longrec = random_seq(200, 3), shortrec = random_seq(10, 4);
    RollingHasher<u128> rh(40, false, DNA, 50);
    rh.for_each_hash([](u128) {}, longrec.data(), longrec.size());
    rh.for_each_hash([](u128) {}, shortrec.data(), shortrec.size());
    REQUIRE(rh.n_in_queue() == 0);
    RollingHasher<uint64_t> rhp(20, false, PROTEIN20, 30);
    rhp.for_each_hash([](uint64_t) {}, longrec.data(), longrec.size());
    rhp.for_each_hash([](uint64_t) {}, shortrec.data(), shortrec.size());
    REQUIRE(rhp.n_in_queue() == 0);
}
