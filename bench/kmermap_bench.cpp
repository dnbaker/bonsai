// CompactKmerMap vs khash_t(c): bytes and classification-shaped lookup
// throughput (reads of K k-mers, a hit/miss mix), printed, never asserted.
// Usage: kmermap_bench [nkmers=50000000] [nreads=200000] [kmers_per_read=150] [hit_frac=0.7] [max_load=0.85]
#include "kmermap.h"

#include <chrono>
#include <cstdio>
#include <random>
#include <vector>

using namespace bns;
using clk = std::chrono::steady_clock;
static double secs(clk::time_point a) {return std::chrono::duration<double>(clk::now() - a).count();}

// khash_t(c) batch lookup with the same prefetch pattern (fair baseline):
// prefetch flags/keys/vals at each key's first probe bucket `ahead` keys early.
template<typename It, typename Ot>
static void kh_get_many_prefetch(const khash_t(c) *h, It f, It l, Ot o, size_t ahead = 16) {
    const khint_t mask = h->n_buckets - 1;
    It pf = f;
    auto pre = [&](u64 k) {
        const khint_t i = kh_int64_hash_func(k) & mask;
        __builtin_prefetch(h->flags + (i >> 4));
        __builtin_prefetch(h->keys + i);
        __builtin_prefetch(h->vals + i);
    };
    for(size_t d = 0; d < ahead && pf != l; ++d) pre(*pf++);
    for(; f != l; ++f, ++o) {
        if(pf != l) pre(*pf++);
        const khint_t ki = kh_get(c, h, *f);
        *o = ki == kh_end(h) ? tax_t(-1): kh_val(h, ki);
    }
}

int main(int argc, char **argv) {
    const size_t n = argc > 1 ? std::strtoull(argv[1], nullptr, 10): 50000000;
    const size_t nreads = argc > 2 ? std::strtoull(argv[2], nullptr, 10): 200000;
    const size_t K = argc > 3 ? std::strtoull(argv[3], nullptr, 10): 150;
    const double hit = argc > 4 ? std::atof(argv[4]): 0.7;
    const double load = argc > 5 ? std::atof(argv[5]): CompactKmerMap::DEFAULT_MAX_LOAD;
    std::mt19937_64 rng(17);
    std::vector<u64> keys(n);
    for(auto &k: keys) k = rng();
    auto tax_of = [](u64 k) {return tax_t((k * 0x9e3779b97f4a7c15ull) >> 40) % 3300000u;};
    auto t = clk::now();
    khash_t(c) *km(kh_init(c));
    int khr;
    for(const u64 k: keys) {const khint_t ki(kh_put(c, km, k, &khr)); kh_val(km, ki) = tax_of(k);}
    const double tkb = secs(t);
    t = clk::now();
    CompactKmerMap cm(kh_size(km), load);
    for(const u64 k: keys) cm.set(k, tax_of(k));
    const double tcb = secs(t);
    std::printf("kmers %zu | khash_t(c) %.0f MB (buckets 2^%d, load %.2f) | CompactKmerMap %.0f MB (q=%u, load %.2f) | ratio %.3f\n",
                cm.size(), khash_bytes(km) / 1e6, __builtin_ctz(km->n_buckets), double(kh_size(km)) / km->n_buckets,
                cm.bytes_used() / 1e6, cm.table().q(), cm.load_factor(), double(cm.bytes_used()) / khash_bytes(km));
    std::printf("build: khash %.2fs, compact %.2fs\n", tkb, tcb);
    // the struct-of-arrays cht_t (4 arrays per hit) for comparison
    sketch::cht::cht_t<uint64_t, sketch::cht::cht_splitmix64_hash, 32> soa(cm.table().q());
    for(const u64 k: keys) soa.assign(k, uint64_t(tax_of(k)) + 1);
    std::printf("cht_t (struct-of-arrays) %.0f MB\n", soa.bytes_used() / 1e6);
    std::vector<u64> q(nreads * K);
    for(auto &x: q) x = (std::uniform_real_distribution<double>(0, 1)(rng) < hit) ? keys[rng() % n]: rng();
    std::vector<tax_t> out(K);
    std::vector<uint64_t> scratch;
    auto run = [&](const char *name, auto &&per_read) {
        u64 sum = 0;
        const auto t0 = clk::now();
        for(size_t r = 0; r < nreads; ++r) {
            per_read(q.data() + r * K, q.data() + (r + 1) * K, out.data());
            for(const tax_t x: out) sum += x;
        }
        const double s = secs(t0);
        std::printf("  %-32s %7.1f M lookups/s  (checksum %llu)\n", name, q.size() / s / 1e6, (unsigned long long)sum);
        return sum;
    };
    std::printf("lookups: %zu reads x %zu k-mers, hit fraction %.2f\n", nreads, K, hit);
    const u64 a = run("khash kh_get (current classify)", [&](const u64 *f, const u64 *l, tax_t *o) {
        for(; f != l; ++f, ++o) {const khint_t ki = kh_get(c, km, *f); *o = ki == kh_end(km) ? tax_t(-1): kh_val(km, ki);}
    });
    const u64 b = run("khash batch + prefetch", [&](const u64 *f, const u64 *l, tax_t *o) {kh_get_many_prefetch(km, f, l, o);});
    const u64 s1 = run("cht_t SoA get_many (prefetch 16)", [&](const u64 *f, const u64 *l, tax_t *o) {
        scratch.resize(l - f); soa.get_many_prefetch(f, l, scratch.begin(), 16);
        for(const uint64_t v: scratch) *o++ = v ? tax_t(v - 1): tax_t(-1);
    });
    const u64 c1 = run("cht_map get (scalar)", [&](const u64 *f, const u64 *l, tax_t *o) {for(; f != l; ++f) *o++ = cm.get(*f);});
    const u64 d = run("cht_map get_many (prefetch 16)", [&](const u64 *f, const u64 *l, tax_t *o) {cm.get_many(f, l, o, scratch, 16);});
    const u64 e = run("cht_map get_many (prefetch 32)", [&](const u64 *f, const u64 *l, tax_t *o) {cm.get_many(f, l, o, scratch, 32);});
    std::printf("checksums %s\n", (a == b && b == s1 && s1 == c1 && c1 == d && d == e) ? "all equal": "DIFFER");
    kh_destroy(c, km);
}
