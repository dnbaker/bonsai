#ifndef BONSAI_KMERMAP_H__
#define BONSAI_KMERMAP_H__
// CompactKmerMap: the k-mer -> taxon (LCA) database in sketch's compact hash
// table (hll/include/sketch/cht.h) instead of khash_t(c).
//
// Layout: cht_map_t<uint64_t, cht_splitmix64_hash, 32> (sketch's single-stream exact map).
// - The 64-bit k-mer is mixed by an invertible SplitMix64. Only 64 - q
//   remainder bits are stored, since the home slot supplies the top q bits.
// - The value is tax + 1 (cht's value 0 means empty). A taxon of tax_t(-1)
//   is not storable.
// - Per slot: ONE packed cell of 64 - q key bits + 32 value bits + 8
//   displacement bits, with value 0 = empty. A probe step is one 16-byte
//   load, so a lookup usually touches a single cache line.
// - khash_t(c) spends 96 bits per bucket plus 2 flag bits, at buckets = the
//   next power of two above n / 0.77.
//
// Lookups are home-exact (get / get_many_prefetch), so there are no merged
// keys and no false positives. Missing k-mers return NOT_FOUND = tax_t(-1),
// the same convention as Database::get_lca.
//
// Classification is lookup-bound: use get_many on all k-mers of a read at
// once. It hashes each k-mer once and prefetches its home lines ahead of use.
// Regression: test/regress_compact_kmermap.cpp. Memory/speed: bench/kmermap_bench.cpp.
#include "util.h"
#include "database.h"
#include "sketch/cht.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>

namespace bns {

class CompactKmerMap {
public:
    using table_type = sketch::cht::cht_map_t<uint64_t, sketch::cht::cht_splitmix64_hash, 32>;
    static constexpr tax_t NOT_FOUND = std::numeric_limits<tax_t>::max();
    static constexpr double DEFAULT_MAX_LOAD = 0.85;

    // Smallest table that holds n entries at load <= max_load (q in [2, 40]).
    static uint32_t q_for(size_t n, double max_load = DEFAULT_MAX_LOAD) {
        uint32_t q = 2;
        while(q < 40 && double(n) > max_load * std::ldexp(1., int(q))) ++q;
        if(double(n) > max_load * std::ldexp(1., int(q)))
            throw std::length_error("CompactKmerMap: more k-mers than a 2^40-slot table holds at this load");
        return q;
    }

    explicit CompactKmerMap(size_t expected_size, double max_load = DEFAULT_MAX_LOAD)
        : t_(q_for(expected_size, max_load)) {}
    explicit CompactKmerMap(table_type &&t): t_(std::move(t)) {}
    explicit CompactKmerMap(const khash_t(c) *map, double max_load = DEFAULT_MAX_LOAD)
        : t_(q_for(kh_size(map), max_load)) {
        for(khint_t ki = kh_begin(map); ki != kh_end(map); ++ki)
            if(kh_exist(map, ki)) set(kh_key(map, ki), kh_val(map, ki));
    }

    // Insert kmer, or overwrite its taxon (kh_put + kh_val semantics).
    void set(u64 kmer, tax_t tax) {
        if(tax == NOT_FOUND)
            throw std::invalid_argument("CompactKmerMap: taxon tax_t(-1) is not storable");
        t_.assign(kmer, uint64_t(tax) + 1);
    }
    tax_t get(u64 kmer) const {
        const uint64_t v = t_.get(kmer);
        return v ? tax_t(v - 1): NOT_FOUND;
    }
    bool contains(u64 kmer) const {return t_.get(kmer) != 0;}
    // Batched lookups in input order (identical to get() on each k-mer).
    // `scratch` avoids a per-call allocation in hot loops.
    template<typename It, typename Ot>
    void get_many(It f, It l, Ot o, std::vector<uint64_t> &scratch, size_t ahead = 16) const {
        scratch.resize(std::distance(f, l));
        t_.get_many_prefetch(f, l, scratch.begin(), ahead);
        for(const uint64_t v: scratch) *o++ = v ? tax_t(v - 1): NOT_FOUND;
    }
    size_t size() const {return t_.size();}
    double load_factor() const {return t_.load_factor();}
    size_t bytes_used() const {return t_.bytes_used();}
    // f(kmer, tax) for every entry, in table order.
    template<typename F>
    void for_each(F &&f) const {
        for(auto it = t_.begin(); it != t_.end(); ++it) f(u64(*it), tax_t(it.mapped() - 1));
    }
    const table_type &table() const {return t_;}
    void write(std::FILE *fp) const {t_.write(fp);}
    static CompactKmerMap read(std::FILE *fp) {return CompactKmerMap(table_type::read(fp));}

private:
    table_type t_;
};

// Compact classification database file: magic "BNSCKM01", then k, w and the
// spacing exactly as Database<khash_t(c)> writes them, then the cht_map_t
// stream. Uncompressed. `bonsai compact` converts a khash database into one.
struct CompactDatabase {
    unsigned k_ = 0, w_ = 0;
    spvec_t s_;
    std::unique_ptr<CompactKmerMap> db_;
    static constexpr char MAGIC[9] = "BNSCKM01";

    static bool is_compact(const char *fn) {
        std::FILE *fp = std::fopen(fn, "rb");
        if(!fp) return false;
        char m[8];
        const bool ret = std::fread(m, 1, 8, fp) == 8 && std::memcmp(m, MAGIC, 8) == 0;
        std::fclose(fp);
        return ret;
    }
    explicit CompactDatabase(const char *fn) {
        std::FILE *fp = std::fopen(fn, "rb");
        if(!fp) LOG_EXIT("Could not open %s for reading.\n", fn);
        char m[8];
        if(std::fread(m, 1, 8, fp) != 8 || std::memcmp(m, MAGIC, 8))
            throw std::runtime_error(std::string("Not a compact bonsai database: ") + fn);
        if(std::fread(&k_, sizeof(k_), 1, fp) != 1 || std::fread(&w_, sizeof(w_), 1, fp) != 1)
            throw std::runtime_error("Error: Could not read k/w from compact database");
        s_ = spvec_t(k_ - 1);
        for(auto &v: s_) {
            uint8_t bv;
            if(std::fread(&bv, sizeof(bv), 1, fp) != 1) throw std::runtime_error("Error: Could not read spacing from file");
            v = bv;
        }
        db_.reset(new CompactKmerMap(CompactKmerMap::read(fp)));
        std::fclose(fp);
    }
    static void write(const char *fn, unsigned k, unsigned w, const spvec_t &s, const CompactKmerMap &m) {
        std::FILE *fp = std::fopen(fn, "wb");
        if(!fp) LOG_EXIT("Could not open %s for writing.\n", fn);
        if(std::fwrite(MAGIC, 1, 8, fp) != 8) throw std::runtime_error("Error writing database");
        if(std::fwrite(&k, sizeof(k), 1, fp) != 1 || std::fwrite(&w, sizeof(w), 1, fp) != 1)
            throw std::runtime_error("Error writing database");
        for(const auto v: s) {
            const uint8_t bv = static_cast<uint8_t>(v);
            if(std::fwrite(&bv, sizeof(bv), 1, fp) != 1) throw std::runtime_error("Error writing database");
        }
        m.write(fp);
        std::fclose(fp);
    }
};

// Bytes khash_t(c) holds for its buckets (keys + vals + flags), for comparison.
inline size_t khash_bytes(const khash_t(c) *map) {
    return size_t(map->n_buckets) * (sizeof(*map->keys) + sizeof(*map->vals))
         + __ac_fsize(map->n_buckets) * sizeof(khint32_t);
}

} // namespace bns

#endif /* BONSAI_KMERMAP_H__ */
