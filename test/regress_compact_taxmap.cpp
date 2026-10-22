// CompactTaxMap (taxmap.h) must be a drop-in for khash_t(p): identical parents
// for every taxid, NOT_FOUND for absent ids, identical node_depth / node_dist /
// lca on random pairs, a complete iteration, and build_compact_parent_map must
// parse files exactly as build_parent_map does (including the compact and
// pipeless rows from regress_parent_map.cpp). Also reports bytes vs khash.
#include "hashmap.h"
#include "test/catch.hpp"

#include "taxmap.h"

#include <cstdio>
#include <random>
#include <string>
#include <vector>

using namespace bns;

namespace {

struct TmpTax {
    std::string path;
    TmpTax(const std::string &content): path([] {
        char tmpl[] = "/tmp/bns_ctax_XXXXXX";
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

// NCBI-like tree: root 1, sparse ids up to ~3e6 plus a few near 2^31, each node's
// parent drawn among earlier nodes with a bias toward recent ones (depth ~10-40).
std::vector<std::pair<tax_t, tax_t>> synth_tree(size_t n, uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::vector<std::pair<tax_t, tax_t>> rows;
    std::vector<tax_t> ids{1};
    bns::uset<tax_t> used{1, 0};
    rows.emplace_back(1, 0);
    while(rows.size() < n) {
        tax_t id = (rng() % 50 == 0) ? tax_t(2000000000u + rng() % 100000000u): tax_t(2 + rng() % 3000000u);
        if(!used.insert(id).second) continue;
        const size_t m = ids.size();
        const size_t pi = (rng() & 1) ? m - 1 - std::min<size_t>(m - 1, rng() % 64): rng() % m;
        rows.emplace_back(id, ids[pi]);
        ids.push_back(id);
    }
    return rows;
}

} // namespace

TEST_CASE("CompactTaxMap matches khash_t(p) on parents, depth, dist and lca") {
    const auto rows = synth_tree(200000, 7);
    khash_t(p) *km(kh_init(p));
    int khr;
    for(const auto &r: rows) {
        const khint_t ki(kh_put(p, km, r.first, &khr)); // separate statement: kh_put may realloc vals
        kh_val(km, ki) = r.second;
    }
    CompactTaxMap cm(km);
    REQUIRE(cm.size() == kh_size(km));
    for(const auto &r: rows) REQUIRE(cm.parent(r.first) == get_parent(km, r.first));
    std::mt19937_64 rng(11);
    size_t absent = 0;
    for(int i = 0; i < 100000; ++i) {
        const tax_t t = tax_t(rng());
        if(kh_get(p, km, t) != kh_end(km)) continue;
        ++absent;
        REQUIRE(cm.parent(t) == CompactTaxMap::NOT_FOUND);
        REQUIRE(!cm.contains(t));
    }
    REQUIRE(absent > 90000);
    for(int i = 0; i < 20000; ++i) {
        const tax_t a = rows[rng() % rows.size()].first, b = rows[rng() % rows.size()].first;
        REQUIRE(lca(cm, a, b) == lca(km, a, b));
        REQUIRE(node_depth(cm, a) == node_depth(km, a));
        const tax_t c = lca(km, a, b);
        if(c != a) REQUIRE(node_dist(cm, a, c) == node_dist(km, a, c));
    }
    std::vector<tax_t> group;
    for(int i = 0; i < 9; ++i) group.push_back(rows[rng() % rows.size()].first);
    REQUIRE(lca(cm, group) == lca(km, group));
    size_t visited = 0;
    cm.for_each([&](tax_t c, tax_t par) {
        ++visited;
        REQUIRE(get_parent(km, c) == par);
    });
    REQUIRE(visited == kh_size(km));
    REQUIRE(cm.table().structure_check());
    std::fprintf(stderr, "CompactTaxMap: %zu taxa, %zu bytes (load %.2f) vs khash_t(p) %zu bytes (%.1f%%)\n",
                 cm.size(), cm.bytes_used(), cm.load_factor(), khash_bytes(km),
                 100. * cm.bytes_used() / khash_bytes(km));
    REQUIRE(cm.bytes_used() < khash_bytes(km));
    kh_destroy(p, km);
}

TEST_CASE("build_compact_parent_map parses exactly like build_parent_map") {
    for(const std::string &content: {std::string("2|131567\n3|\t2|\n131567 | 1\n"),
                                    std::string("2\t|\t1\ngarbage line with no pipe\n3\t|\t2\n"),
                                    std::string("# comment\n\n5\t|\t1\t|\tgenus\n5\t|\t2\n2\t|\t1\n1\t|\t1\n")}) {
        TmpTax tax(content);
        khash_t(p) *km(build_parent_map(tax.path.data()));
        CompactTaxMap cm(build_compact_parent_map(tax.path.data()));
        REQUIRE(cm.size() == kh_size(km));
        for(khint_t ki = kh_begin(km); ki != kh_end(km); ++ki)
            if(kh_exist(km, ki)) REQUIRE(cm.parent(kh_key(km, ki)) == kh_val(km, ki));
        REQUIRE(cm.parent(1) == 0);                  // root forced to 0
        REQUIRE(cm.parent(0) == CompactTaxMap::NOT_FOUND);
        kh_destroy(p, km);
    }
}

TEST_CASE("TaxIdRevHash is a bijection (hash o invert == identity)") {
    TaxIdRevHash h;
    std::mt19937_64 rng(3);
    for(uint64_t x: {uint64_t(0), uint64_t(1), uint64_t(0xffffffffu), uint64_t(0x80000000u)}) {
        REQUIRE(h.invert(h.hash(x)) == x);
        REQUIRE(h.hash(h.invert(x)) == x);
    }
    for(int i = 0; i < 1000000; ++i) {
        const uint64_t x = uint32_t(rng());
        REQUIRE(h.invert(h.hash(x)) == x);
    }
}

TEST_CASE("CompactTaxMap rejects tax_t(-1) parents and overwrites like kh_put") {
    CompactTaxMap cm(16);
    cm.set_parent(10, 1);
    cm.set_parent(10, 2);
    REQUIRE(cm.parent(10) == 2);
    REQUIRE(cm.size() == 1);
    REQUIRE_THROWS_AS(cm.set_parent(11, tax_t(-1)), std::invalid_argument);
}
