#ifndef BONSAI_KMERMAP_H__
#define BONSAI_KMERMAP_H__
// CompactKmerMap: the k-mer -> taxon (LCA) database in sketch's compact hash
// table (hll/include/sketch/cht.h) instead of khash_t(c).
//
// Layout: cht_map_t<uint64_t, cht_splitmix64_hash, 0> (sketch's single-stream
// exact map, value width chosen at construction).
// - The 64-bit k-mer is mixed by an invertible SplitMix64. Only 64 - q
//   remainder bits are stored, since the home slot supplies the top q bits.
// - DENSE TAXON IDS: a database holds far fewer distinct taxa than 2^32, so
//   each value is a dense id + 1 (0 means empty). A small table maps dense
//   ids back to taxids (taxa_; ids follow sorted taxid order, so one
//   database always encodes the same way). The value width is the minimum
//   for the taxon count: ceil(log2(ntaxa + 1)) bits, e.g. 20 for ~1M taxa,
//   instead of 32. The taxon table is a few MB and stays in cache.
// - Per slot: ONE packed cell of 64 - q key bits + value bits + 8
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
#include <algorithm>
#include <memory>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace bns {

class CompactKmerMap {
public:
    using table_type = sketch::cht::cht_map_t<uint64_t, sketch::cht::cht_splitmix64_hash, 0>;
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
    // Value width for dense ids 1..ntaxa.
    static uint32_t bits_for(size_t ntaxa) {
        uint32_t b = 1;
        while(b < 48 && ((uint64_t(1) << b) - 1) < uint64_t(ntaxa)) ++b;
        return b;
    }

    // Incremental build: room for expected_size k-mers and up to max_taxa
    // distinct taxa (set() throws past either).
    CompactKmerMap(size_t expected_size, size_t max_taxa, double max_load = DEFAULT_MAX_LOAD)
        : t_(q_for(expected_size, max_load), bits_for(max_taxa)), max_taxa_(max_taxa) {}
    // Conversion from khash_t(c): one pass to collect the distinct taxa (the
    // minimum value width), one to fill.
    explicit CompactKmerMap(const khash_t(c) *map, double max_load = DEFAULT_MAX_LOAD)
        : CompactKmerMap(map, collect_taxa(map), max_load) {}

    // Insert kmer, or overwrite its taxon (kh_put + kh_val semantics).
    void set(u64 kmer, tax_t tax) {
        if(tax == NOT_FOUND)
            throw std::invalid_argument("CompactKmerMap: taxon tax_t(-1) is not storable");
        t_.assign(kmer, dense_id(tax));
    }
    tax_t get(u64 kmer) const {
        const uint64_t v = t_.get(kmer);
        return v ? taxa_[v - 1]: NOT_FOUND;
    }
    bool contains(u64 kmer) const {return t_.get(kmer) != 0;}
    // Batched lookups in input order (identical to get() on each k-mer).
    // `scratch` avoids a per-call allocation in hot loops.
    template<typename It, typename Ot>
    void get_many(It f, It l, Ot o, std::vector<uint64_t> &scratch, size_t ahead = 16) const {
        scratch.resize(std::distance(f, l));
        t_.get_many_prefetch(f, l, scratch.begin(), ahead);
        for(const uint64_t v: scratch) *o++ = v ? taxa_[v - 1]: NOT_FOUND;
    }
    size_t size() const {return t_.size();}
    size_t num_taxa() const {return taxa_.size();}
    uint32_t value_bits() const {return t_.vbits();}
    double load_factor() const {return t_.load_factor();}
    // table + taxon table (+ the build-time id index while it exists)
    size_t bytes_used() const {
        return t_.bytes_used() + taxa_.size() * sizeof(tax_t)
             + ids_.size() * (sizeof(tax_t) + sizeof(uint32_t) + 2 * sizeof(void *));
    }
    // Drop the taxid -> dense id index (needed only by set()); it is rebuilt
    // on demand if set() is called again.
    void finalize() {decltype(ids_)().swap(ids_);}
    // f(kmer, tax) for every entry, in table order.
    template<typename F>
    void for_each(F &&f) const {
        for(auto it = t_.begin(); it != t_.end(); ++it) f(u64(*it), taxa_[it.mapped() - 1]);
    }
    const table_type &table() const {return t_;}
    void write(std::FILE *fp) const {
        const uint64_t nt = taxa_.size(), mt = max_taxa_;
        if(std::fwrite(&nt, sizeof(nt), 1, fp) != 1 || std::fwrite(&mt, sizeof(mt), 1, fp) != 1
           || (nt && std::fwrite(taxa_.data(), sizeof(tax_t), nt, fp) != nt))
            throw std::runtime_error("CompactKmerMap: write failed");
        t_.write(fp);
    }
    static CompactKmerMap read(std::FILE *fp) {
        uint64_t nt, mt;
        if(std::fread(&nt, sizeof(nt), 1, fp) != 1 || std::fread(&mt, sizeof(mt), 1, fp) != 1)
            throw std::runtime_error("CompactKmerMap: truncated taxon table header");
        std::vector<tax_t> taxa(nt);
        if(nt && std::fread(taxa.data(), sizeof(tax_t), nt, fp) != nt)
            throw std::runtime_error("CompactKmerMap: truncated taxon table");
        CompactKmerMap ret(table_type::read(fp));
        if(ret.t_.vbits() < bits_for(nt))
            throw std::runtime_error("CompactKmerMap: value width too small for the taxon table");
        ret.taxa_ = std::move(taxa);
        ret.max_taxa_ = mt;
        return ret;
    }

private:
    table_type t_;
    std::vector<tax_t> taxa_;                       // dense id - 1 -> taxid
    std::unordered_map<tax_t, uint32_t> ids_;       // taxid -> dense id (build-time index)
    size_t max_taxa_ = 0;

    explicit CompactKmerMap(table_type &&t): t_(std::move(t)) {}

    // Distinct taxa of a khash_t(c) with a taxid -> dense id lookup for the
    // fill pass. Taxids below DIRECT_MAX (NCBI ids reach ~3.3M) use a direct
    // array, so the per-k-mer lookup is one array read; larger custom ids
    // fall back to a hash map. Dense ids follow sorted taxid order either way.
    struct TaxIndex {
        static constexpr tax_t DIRECT_MAX = tax_t(1) << 26;   // <= 256 MiB of u32 at the cap
        std::vector<tax_t> taxa;                             // sorted distinct taxa
        std::vector<uint32_t> direct;                        // taxid -> dense id (when max < DIRECT_MAX)
        std::unordered_map<tax_t, uint32_t> hashed;          // fallback
        uint32_t id(tax_t t) const {return direct.empty() ? hashed.at(t): direct[t];}
    };
    static TaxIndex collect_taxa(const khash_t(c) *map) {
        TaxIndex ix;
        tax_t mx = 0;
        for(khint_t ki = kh_begin(map); ki != kh_end(map); ++ki)
            if(kh_exist(map, ki)) mx = std::max(mx, kh_val(map, ki));
        if(mx == NOT_FOUND) throw std::invalid_argument("CompactKmerMap: taxon tax_t(-1) is not storable");
        if(mx < TaxIndex::DIRECT_MAX) {
            ix.direct.assign(size_t(mx) + 1, 0);
            for(khint_t ki = kh_begin(map); ki != kh_end(map); ++ki)
                if(kh_exist(map, ki)) ix.direct[kh_val(map, ki)] = 1;
            uint32_t next = 0;
            for(size_t t = 0; t <= mx; ++t)
                if(ix.direct[t]) {ix.direct[t] = ++next; ix.taxa.push_back(tax_t(t));}
        } else {
            std::unordered_set<tax_t> seen;
            for(khint_t ki = kh_begin(map); ki != kh_end(map); ++ki)
                if(kh_exist(map, ki) && seen.insert(kh_val(map, ki)).second) ix.taxa.push_back(kh_val(map, ki));
            std::sort(ix.taxa.begin(), ix.taxa.end());
            ix.hashed.reserve(ix.taxa.size());
            for(size_t i = 0; i < ix.taxa.size(); ++i) ix.hashed.emplace(ix.taxa[i], uint32_t(i + 1));
        }
        return ix;
    }
    CompactKmerMap(const khash_t(c) *map, TaxIndex &&ix, double max_load)
        : t_(q_for(kh_size(map), max_load), bits_for(ix.taxa.size())) {
        for(khint_t ki = kh_begin(map); ki != kh_end(map); ++ki)
            if(kh_exist(map, ki)) t_.assign(kh_key(map, ki), ix.id(kh_val(map, ki)));
        taxa_ = std::move(ix.taxa);
        max_taxa_ = taxa_.size();
    }
    uint32_t dense_id(tax_t tax) {
        if(ids_.empty() && !taxa_.empty())
            for(size_t i = 0; i < taxa_.size(); ++i) ids_.emplace(taxa_[i], uint32_t(i + 1));
        const auto it = ids_.find(tax);
        if(it != ids_.end()) return it->second;
        if(taxa_.size() >= max_taxa_ || taxa_.size() + 1 > t_.sat())
            throw std::length_error("CompactKmerMap: more distinct taxa than this map was sized for");
        taxa_.push_back(tax);
        return ids_.emplace(tax, uint32_t(taxa_.size())).first->second;
    }
};

// Compact classification database file: magic "BNSCKM02", then k, w and the
// spacing exactly as Database<khash_t(c)> writes them, then the dense taxon
// table (count, capacity, taxids) and the cht_map_t stream. Uncompressed. `bonsai compact` converts a khash database into one.
struct CompactDatabase {
    unsigned k_ = 0, w_ = 0;
    spvec_t s_;
    std::unique_ptr<CompactKmerMap> db_;
    static constexpr char MAGIC[9] = "BNSCKM02";

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
