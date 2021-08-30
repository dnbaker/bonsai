#include "test/catch.hpp"
#include "util.h"
using namespace bns;

#define is_pow2(x) ((x & (x - 1)) == 0)

TEST_CASE("KhashSerial") {
    khash_t(c) *th(kh_init(c)), *ti(nullptr);
    khint_t ki;
    int khr;
    for(size_t i(0); i < 128; ++i) {
        uint64_t k((i << 14) | (i + 2));
        int val(i);
        ki = kh_put(c, th, k, &khr);
        kh_val(th, ki) = val;
    }
    for(khiter_t ki = 0; ki < th->n_buckets; ++ki) {
        if(kh_exist(th, ki)) {
            std::fprintf(stderr, "Key %zu->%u\n", size_t(th->keys[ki]), int(th->vals[ki]));
        }
    }
    khash_write(th, "__zomg__");
    ti = khash_load<khash_t(c)>("__zomg__");
    REQUIRE(system("rm __zomg__") == 0);
    for(size_t i(0); i < th->n_buckets; ++i) {
        REQUIRE(th->flags[__ac_fsize(i)] == ti->flags[__ac_fsize(i)]);
        REQUIRE(th->vals[i] == ti->vals[i]);
        REQUIRE(th->keys[i] == ti->keys[i]);
    }
    //std::fclose(fp);
    kh_destroy(c, th);
    kh_destroy(c, ti);
}

TEST_CASE("roundup64") {
    for(size_t i(0); i < 1 << 10; ++i) {
        size_t d(((uint64_t)rand() << 32) | rand());
        if(is_pow2(d)) --d;
        REQUIRE(__builtin_clzll(d) - 1 == __builtin_clzll(roundup64(d)));
        if(is_pow2(d)) continue;
        REQUIRE(__builtin_clzll(d) - 1 == __builtin_clzll(roundup64(d)));
    }
}
