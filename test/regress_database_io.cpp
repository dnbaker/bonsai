// Database binary IO: the plain and gzipped write/read paths must round-trip
// k, w, the spacing vector, and every key/value pair. (The spacing return-value
// checks once compared fwrite/fread item counts to byte counts, throwing for
// k != 2, and the gzipped reader mixed stdio with raw read(2), which drops
// prefetched bytes on pipes.)
#include "test/catch.hpp"

#include "database.h"
#include "util.h"

#include <cstdio>
#include <string>
#include <vector>

using namespace bns;

namespace {

struct TmpDb {
    std::string path;
    explicit TmpDb(const std::string &suffix): path([] {
        char tmpl[] = "/tmp/bns_db_XXXXXX";
        const char *ret(mkdtemp(tmpl));
        REQUIRE(ret != nullptr);
        return std::string(ret);
    }() + suffix) {}
    ~TmpDb() { std::remove(path.data()); }
};

void fill_db(Database<khash_t(c)> &db) {
    kh_resize(c, db.db_, 64);
    for(u64 k(1); k < 64; ++k) {
        int khr;
        auto it = kh_put(c, db.db_, k * 2654435761u, &khr);
        REQUIRE(khr >= 0);
        kh_val(db.db_, it) = tax_t(k * 3);
    }
}

void check_db(const Database<khash_t(c)> &db) {
    REQUIRE(db.k_ == 11);
    REQUIRE(db.w_ == 11);
    REQUIRE(db.s_.size() == 10);
    for(const auto v: db.s_) REQUIRE(v == 0); // unspaced differences
    REQUIRE(kh_size(db.db_) == 63);
    for(u64 k(1); k < 64; ++k) {
        const auto ki = kh_get(c, db.db_, k * 2654435761u);
        REQUIRE(ki != kh_end(db.db_));
        REQUIRE(kh_val(db.db_, ki) == tax_t(k * 3));
    }
    // get_lca reports the value for present keys and -1u for absent ones.
    REQUIRE(db.get_lca(2654435761u) == 3);
    REQUIRE(db.get_lca(42) == -1u);
}

} // namespace

TEST_CASE("Database binary IO round-trips plain and gzipped") {
    const Spacer sp(11, 11);
    SECTION("plain") {
        TmpDb tmp(".db");
        {
            Database<khash_t(c)> db(sp, 1, kh_init(c));
            fill_db(db);
            db.write(tmp.path.data(), false);
        }
        Database<khash_t(c)> db(tmp.path.data());
        check_db(db);
    }
    SECTION("gzipped") {
        TmpDb tmp(".db.gz");
        {
            Database<khash_t(c)> db(sp, 1, kh_init(c));
            fill_db(db);
            db.write(tmp.path.data(), true);
        }
        Database<khash_t(c)> db(tmp.path.data());
        check_db(db);
    }
    SECTION("spaced") {
        TmpDb tmp(".db");
        const spvec_t spaces {0, 1, 0, 2, 0, 1, 0, 2, 0, 1};
        const Spacer ssp(11, 11, spaces);
        {
            Database<khash_t(c)> db(ssp, 1, kh_init(c));
            fill_db(db);
            db.write(tmp.path.data(), false);
        }
        Database<khash_t(c)> db(tmp.path.data());
        REQUIRE(db.k_ == 11);
        REQUIRE(db.s_.size() == 10);
        REQUIRE(db.s_ == spaces);
        REQUIRE(kh_size(db.db_) == 63);
    }
}
