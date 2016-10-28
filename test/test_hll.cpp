#include "test/catch.hpp"
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <string>
#include "lib/hll.h"
#include "lib/encoder.h"

using namespace kpg;

int within_bounds(hll_t &hll, double val) {
    return std::abs(hll.report() - val) <= hll.est_err();
}

std::vector<std::pair<uint64_t, uint64_t>> pairs {
    {8, 16},
    {16, 27},
    //{20, 32},
    //{23, 36},
    //{27, 40},
};
std::vector<std::string> paths {
    "test/GCF_000302455.1_ASM30245v1_genomic.fna.gz",
    "test/GCF_000762265.1_ASM76226v1_genomic.fna.gz",
    "test/GCF_000953115.1_DSM1535_genomic.fna.gz",
    "test/GCF_001723155.1_ASM172315v1_genomic.fna.gz"
};

TEST_CASE("hll") {
    //for(auto x: {8, 12, 16, 20, 24, 28, 32, 36, 42}) {
    for(auto &pair: pairs) {
        size_t val(1ull << pair.second);
        size_t hsz(pair.first);
        hll_t hll(hsz);
        do {
            if(hll.is_ready() && hll.m_ != hsz) {
                --hsz;
                size_t inc((int)(std::log2(val / hll.report()) + 0.5) / 2 + 1);
                if(!inc) ++inc;
                hsz += inc;
                fprintf(stderr, "Estimated bound surrounding %lf of %lf is less than found %lf. Incrementing hsz by %zu to %zu\n",
                        hll.report(), hll.est_err(), std::abs(hll.report() - (double)val), inc, hsz);
                hll = hll_t(hsz);
            }
            for(size_t i(0); i < val; ++i) {
                hll.add(u64hash(((uint64_t)rand() << 32) | rand()));
            }
        } while(0);//while(!within_bounds(hll, val) && ++hsz);
        //fprintf(stderr, "for x = %zu and value is %zu, est is %lf, with %lf. Within bounds for hsz %zu? %s\n",
        //        pair.second, val, hll.report(), hll.est_err(), hsz, within_bounds(hll, val) ? "true": "false");
        REQUIRE(within_bounds(hll, val));
    }
}

TEST_CASE("parallel hll vs hash") {
    spvec_t vec;
    static const size_t np(25);
    ssize_t exact = count_cardinality<lex_score>(paths, 31, 31, vec, nullptr, 2);
    ssize_t inexact = estimate_cardinality<lex_score, np>(paths, 31, 31, vec, nullptr, 2);
    REQUIRE(std::abs(exact - inexact) / exact < 1.03896 * exact / std::sqrt(1 << np));
}
