#include "test/catch.hpp"
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <string>
#include <unordered_set>
#include "hll/hll.h"
#include "lib/encoder.h"

using namespace emp;

// This file tests pure hll functionality
// and cardinal estimation, including
// a test for multi-threaded use.

int within_bounds(hll::hll_t &hll, double val) {
    return std::abs(hll.report() - val) <= hll.est_err();
}

std::vector<std::string> paths {
    "test/GCF_000302455.1_ASM30245v1_genomic.fna.gz",
    "test/GCF_000762265.1_ASM76226v1_genomic.fna.gz",
    "test/GCF_000953115.1_DSM1535_genomic.fna.gz",
    "test/GCF_001723155.1_ASM172315v1_genomic.fna.gz"
};

TEST_CASE("hll") {
    std::vector<std::pair<uint64_t, uint64_t>> pairs {
        //{24, 40},
        {22, 26},
        {19, 24},
        {15, 21}
    };
    //for(auto x: {8, 12, 16, 20, 24, 28, 32, 36, 42}) {
    for(const auto &pair: pairs) {
        const size_t val(1ull << pair.second), hsz(pair.first);
        hll::hll_t hll(hsz);
        for(size_t i(0); i < val; ++i) {
            hll.add(wang_hash(((uint64_t)rand() << 32) | rand()));
        }
        hll.sum();
        if(std::abs(hll.report() - val) > hll.est_err()) {
            fprintf(stderr, "Warning: Failing for %u, %u.\n", (unsigned)pair.first, (unsigned)pair.second);
        }
        REQUIRE(std::abs(hll.report() - val) <= hll.est_err());
    }
}

TEST_CASE("phll") {
    spvec_t vec;
    //static const size_t nps[] {23, 24, 25, 29, 30, 31};
    static const size_t nps[] {23};
    for(const auto np: nps) {
        const ssize_t exact(count_cardinality<lex_score>(paths, 31, 31, vec, nullptr, 2));
        const ssize_t inexact(estimate_cardinality<lex_score>(paths, 31, 31, vec, nullptr, 2, np));
        fprintf(stderr, "For np %zu, we have %lf for expected and %lf for measured as correct as we expect to be, theoretically for %zu and %zu.\n",
                np, (1.03896 / std::sqrt(1ull << np)), (std::abs((double)exact - inexact) / inexact), exact, inexact);
        REQUIRE(std::abs((double)exact - inexact) < 1.03896 * inexact / std::sqrt(1uLL << np));
    }
}

TEST_CASE("Invert Wang") {
    for(size_t i{0}; i < 10000uL; ++i) {
        REQUIRE(irving_inv_hash(wang_hash(i)) == i);
        uint64_t r(rand64());
        REQUIRE(irving_inv_hash(wang_hash(r)) == r);
    }
}
