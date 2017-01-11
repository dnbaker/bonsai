#include "test/catch.hpp"
#include "lib/database.h"
#include <stdlib.h>
using namespace kpg;

#if 1

TEST_CASE("STUFF") {
    khash_t(c) *h(kh_init(c));
    khiter_t ki;
    int khr;
    for(size_t i(0); i < 100000; ++i) {
        ki = kh_put(c, h, rand() * i, &khr);
        kh_val(h, ki) = rand();
    }
    spvec_t v(30, 1);
    Database<khash_t(c)> db(31u, 60u, v, 1, h);
    fprintf(stderr, "Made db!\n");
    db.write("test.db");
    fprintf(stderr, "Written!\n");
    Database<khash_t(c)> read("test.db");
    LOG_DEBUG("Read the stuff from disk\n");
    REQUIRE(read.k_ == db.k_);
    REQUIRE(read.w_ == db.w_);
    for(unsigned i(0); i < db.s_.size(); ++i) {
        REQUIRE(db.s_[i] == read.s_[i]);
    }
    for(unsigned i(0); i != kh_end(read.db_); ++i) {
        REQUIRE(kh_key(read.db_, i) == kh_key(h, i));
        REQUIRE(kh_val(read.db_, i) == kh_val(h, i));
        REQUIRE(kh_exist(read.db_, i) == kh_exist(h, i));
    }
}

#endif
