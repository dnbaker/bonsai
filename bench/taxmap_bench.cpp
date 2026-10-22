// CompactTaxMap vs khash_t(p): bytes, parent-lookup and lca throughput on an
// NCBI-scale synthetic taxonomy (printed, never asserted).
// Usage: taxmap_bench [ntaxa=2700000] [max_load=0.85]
#include "hashmap.h"
#include "taxmap.h"

#include <chrono>
#include <cstdio>
#include <random>
#include <unordered_set>
#include <vector>

using namespace bns;
using clk = std::chrono::steady_clock;
static double secs(clk::time_point a) {return std::chrono::duration<double>(clk::now() - a).count();}

int main(int argc, char **argv) {
    const size_t n = argc > 1 ? std::strtoull(argv[1], nullptr, 10): 2700000;
    const double load = argc > 2 ? std::atof(argv[2]): CompactTaxMap::DEFAULT_MAX_LOAD;
    std::mt19937_64 rng(13);
    std::vector<tax_t> ids{1};
    std::vector<std::pair<tax_t, tax_t>> rows{{1, 0}};
    bns::uset<tax_t> used{0, 1};
    while(rows.size() < n) {   // NCBI-like: ids up to ~3.3e6, a few large custom ids
        const tax_t id = (rng() % 100 == 0) ? tax_t(2000000000u + rng() % 100000000u): tax_t(2 + rng() % 3300000u);
        if(!used.insert(id).second) continue;
        const size_t m = ids.size();
        rows.emplace_back(id, ids[(rng() & 1) ? m - 1 - std::min<size_t>(m - 1, rng() % 64): rng() % m]);
        ids.push_back(id);
    }
    auto t = clk::now();
    khash_t(p) *km(kh_init(p));
    int khr;
    for(const auto &r: rows) {const khint_t ki(kh_put(p, km, r.first, &khr)); kh_val(km, ki) = r.second;}
    const double tk_build = secs(t);
    t = clk::now();
    CompactTaxMap cm(rows.size(), load);
    for(const auto &r: rows) cm.set_parent(r.first, r.second);
    const double tc_build = secs(t);
    std::printf("taxa %zu | khash_t(p) %.1f MB (buckets %u) | CompactTaxMap %.1f MB (load %.2f, q=%u) | ratio %.3f\n",
                cm.size(), khash_bytes(km) / 1e6, km->n_buckets, cm.bytes_used() / 1e6, cm.load_factor(),
                cm.table().q(), double(cm.bytes_used()) / khash_bytes(km));
    std::printf("build: khash %.3fs, compact %.3fs\n", tk_build, tc_build);
    std::vector<tax_t> q(4000000);
    for(auto &x: q) x = (rng() & 3) ? rows[rng() % rows.size()].first: tax_t(rng());   // 75% hits
    uint64_t sk = 0, sc = 0;
    t = clk::now(); for(const tax_t x: q) sk += get_parent(km, x); const double tk = secs(t);
    t = clk::now(); for(const tax_t x: q) sc += get_parent(cm, x); const double tc = secs(t);
    std::printf("get_parent x%zu: khash %.1f M/s, compact %.1f M/s (checksums %s)\n", q.size(),
                q.size() / tk / 1e6, q.size() / tc / 1e6, sk == sc ? "equal": "DIFFER");
    std::vector<std::pair<tax_t, tax_t>> pairs(300000);
    for(auto &pr: pairs) pr = {rows[rng() % rows.size()].first, rows[rng() % rows.size()].first};
    sk = sc = 0;
    t = clk::now(); for(const auto &pr: pairs) sk += lca(km, pr.first, pr.second); const double lk = secs(t);
    t = clk::now(); for(const auto &pr: pairs) sc += lca(cm, pr.first, pr.second); const double lc = secs(t);
    std::printf("lca x%zu: khash %.2f M/s, compact %.2f M/s (checksums %s)\n", pairs.size(),
                pairs.size() / lk / 1e6, pairs.size() / lc / 1e6, sk == sc ? "equal": "DIFFER");
    kh_destroy(p, km);
}
