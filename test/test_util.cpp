#include "test/catch.hpp"
#include "lib/util.h"
using namespace kpg;

#define is_pow2(x) ((x & (x - 1)) == 0)

TEST_CASE("Khash writes and reads correcly") {
    khash_t(c) *th(kh_init(c)), *ti((khash_t(c) *)calloc(1, sizeof(khash_t(c))));
    khint_t ki;
    int khr;
    for(size_t i(0); i < 128; ++i) {
        uint64_t k((i << 14) | (i + 2));
        int val(i);
        ki = kh_put(c, th, k, &khr);
        kh_val(th, ki) = val;
    }
    write_khash_map(th, "zomg");
    ti = load_khash_map<khash_t(c)>("zomg");
    system("rm zomg");
    for(size_t i(0); i < th->n_buckets; ++i) {
        REQUIRE(th->flags[__ac_fsize(i)] == ti->flags[__ac_fsize(i)]);
        REQUIRE(th->vals[i] == ti->vals[i]);
        REQUIRE(th->keys[i] == ti->keys[i]);
    }
    //fclose(fp);
    kh_destroy(c, th);
    kh_destroy(c, ti);
}

TEST_CASE("roundup64") {
    for(size_t i(0); i < 1 << 4; ++i) {
        size_t d(((uint64_t)rand() << 32) | rand());
        if(is_pow2(d)) d = i;
        REQUIRE(__builtin_clzll(d) - 1 == __builtin_clzll(roundup64(d)));
        if(is_pow2(d)) continue;
        REQUIRE(__builtin_clzll(d) - 1 == __builtin_clzll(roundup64(d)));
    }
}
