// build_parent_map row parsing: compact "child|parent" rows must not be
// silently mangled (the old atoi(pipe + 2) read "100|123" as parent 23 and
// made "100|1" a spurious second root), and pipeless rows must be skipped,
// not inserted with parent tax_t(-1).
#include "test/catch.hpp"

#include "util.h"

#include <cstdio>
#include <string>

using namespace bns;

namespace {

struct TmpTax {
    std::string path;
    TmpTax(const std::string &content): path([] {
        char tmpl[] = "/tmp/bns_tax_XXXXXX";
        const char *ret(mkdtemp(tmpl));
        REQUIRE(ret != nullptr);
        return std::string(ret);
    }() + std::string("/tax.txt")) {
        std::FILE *fp(std::fopen(path.data(), "w"));
        REQUIRE(fp != nullptr);
        std::fputs(content.data(), fp);
        std::fclose(fp);
    }
    ~TmpTax() { std::remove(path.data()); }
};

} // namespace

TEST_CASE("build_parent_map parses compact and nodes.dmp-style rows") {
    TmpTax tax("2|131567\n3|\t2|\n131567 | 1\n");
    khash_t(p) *map(build_parent_map(tax.path.data()));
    REQUIRE(kh_get(p, map, 2) != kh_end(map));
    REQUIRE(kh_val(map, kh_get(p, map, 2)) == 131567);
    REQUIRE(kh_val(map, kh_get(p, map, 3)) == 2);
    REQUIRE(kh_val(map, kh_get(p, map, 131567)) == 1);
    // Root of the tree is forced to 0.
    REQUIRE(kh_val(map, kh_get(p, map, 1)) == 0);
    kh_destroy(p, map);
}

TEST_CASE("build_parent_map skips pipeless rows instead of inserting tax_t(-1) parents") {
    TmpTax tax("2\t|\t1\ngarbage line with no pipe\n3\t|\t2\n");
    khash_t(p) *map(build_parent_map(tax.path.data()));
    REQUIRE(kh_get(p, map, 3) != kh_end(map));
    REQUIRE(kh_val(map, kh_get(p, map, 3)) == 2);
    const khint_t ki(kh_get(p, map, 0));
    REQUIRE(ki == kh_end(map)); // pipeless row never entered the map
    for(khint_t it(kh_begin(map)); it != kh_end(map); ++it)
        if(kh_exist(map, it)) REQUIRE(kh_val(map, it) != tax_t(-1));
    kh_destroy(p, map);
}
