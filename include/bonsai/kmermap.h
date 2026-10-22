#ifndef BONSAI_KMERMAP_H__
#define BONSAI_KMERMAP_H__
// CompactKmerMap: the k-mer -> taxon (LCA) database in sketch's compact hash
// table (hll/include/sketch/cht.h) instead of khash_t(c).
//
// Layout: cht_map_t<uint64_t, cht_mixbits_hash, 0> (sketch's single-stream
// exact map; key and value widths chosen at construction).
// - KEY WIDTH: a 2-bit-encoded k-mer spans 2k bits, so the key is mixed by
//   an invertible hash on exactly that domain (cht_mixbits_hash(2k), k <= 32)
//   and only 2k - q remainder bits are stored, since the home slot supplies
//   the top q bits. With the default of 64 any 64-bit key is accepted. Keys
//   outside the domain are rejected by set() and reported absent by
//   get/get_many: the hash masks its input, so they must not alias.
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
#include "hashmap.h"
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

namespace kmermap_detail {
template<typename KeyT> struct KeyHasher;
template<> struct KeyHasher<u64> {using type = sketch::cht::cht_mixbits_hash; static constexpr unsigned max_bits = 64;};
template<> struct KeyHasher<u128> {using type = sketch::cht::cht_mixbits128_hash; static constexpr unsigned max_bits = 128;};
} // namespace kmermap_detail

// Shard of a key for sharded databases (build -x P): a splitmix64 finalizer,
// independent of the table's own hash, then the high bits times P (uniform).
inline unsigned kmer_shard(u64 k, unsigned nshards) {
    k ^= k >> 30; k *= 0xbf58476d1ce4e5b9ULL;
    k ^= k >> 27; k *= 0x94d049bb133111ebULL;
    k ^= k >> 31;
    return unsigned((static_cast<u128>(k) * nshards) >> 64);
}
inline unsigned kmer_shard(u128 k, unsigned nshards) {
    return kmer_shard(u64(k) ^ (u64(k >> 64) * 0x9e3779b97f4a7c15ULL), nshards);
}

// KeyT = u64: k-mers of up to 64 bits (DNA k <= 32) or 64-bit hashes;
// KeyT = u128: 65..128-bit keys (DNA 32 < k <= 64 exactly, or wide hashes).
// CompactKmerMap / CompactKmerMap128 below.
template<typename KeyT>
class BasicCompactKmerMap {
public:
    using key_type = KeyT;
    using hasher_type = typename kmermap_detail::KeyHasher<KeyT>::type;
    static constexpr unsigned MAX_KEY_BITS = kmermap_detail::KeyHasher<KeyT>::max_bits;
    using table_type = sketch::cht::cht_map_t<KeyT, hasher_type, 0>;
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
    BasicCompactKmerMap(size_t expected_size, size_t max_taxa, double max_load = DEFAULT_MAX_LOAD, unsigned key_bits = MAX_KEY_BITS)
        : t_(q_for(expected_size, max_load), bits_for(max_taxa), hasher_type(key_bits)), max_taxa_(max_taxa) {}
    // Conversion from khash_t(c): one pass to collect the distinct taxa (the
    // minimum value width), one to fill.
    explicit BasicCompactKmerMap(const khash_t(c) *map, double max_load = DEFAULT_MAX_LOAD, unsigned key_bits = 64)
        : BasicCompactKmerMap(map, collect_taxa(map), max_load, key_bits) {}
    // Bits the keys of a khash_t(c) actually span (an upper bound for key_bits).
    static unsigned key_bits_needed(const khash_t(c) *map) {
        u64 acc = 0;
        for(khint_t ki = kh_begin(map); ki != kh_end(map); ++ki) if(kh_exist(map, ki)) acc |= kh_key(map, ki);
        unsigned b = 4;
        while(b < 64 && (acc >> b)) ++b;
        return b;
    }

    // Insert kmer, or overwrite its taxon (kh_put + kh_val semantics).
    void set(KeyT kmer, tax_t tax) {
        if(tax == NOT_FOUND)
            throw std::invalid_argument("CompactKmerMap: taxon tax_t(-1) is not storable");
        if(!in_domain(kmer))
            throw std::invalid_argument("CompactKmerMap: k-mer outside the key domain (key_bits)");
        t_.assign(kmer, dense_id(tax));
    }
    tax_t get(KeyT kmer) const {
        if(!in_domain(kmer)) return NOT_FOUND;
        const uint64_t v = t_.get(kmer);
        return v ? taxa_[v - 1]: NOT_FOUND;
    }
    bool contains(KeyT kmer) const {return in_domain(kmer) && t_.get(kmer) != 0;}
    bool in_domain(KeyT kmer) const {return (kmer & ~mask_) == 0;}
    unsigned key_bits() const {return t_.kbits();}
    // Batched lookups in input order (identical to get() on each k-mer).
    // `scratch` avoids a per-call allocation in hot loops.
    template<typename It, typename Ot>
    void get_many(It f, It l, Ot o, std::vector<uint64_t> &scratch, size_t ahead = 16) const {
        scratch.resize(std::distance(f, l));
        t_.get_many_prefetch(f, l, scratch.begin(), ahead);
        if(mask_ == ~KeyT(0)) {
            for(const uint64_t v: scratch) *o++ = v ? taxa_[v - 1]: NOT_FOUND;
        } else {
            size_t i = 0;
            for(It it = f; it != l; ++it, ++i)
                *o++ = (scratch[i] && in_domain(*it)) ? taxa_[scratch[i] - 1]: NOT_FOUND;
        }
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
    // Pre-assign dense ids 1..n to `taxa` in the given order (e.g. sorted, so
    // a database encodes deterministically); only on an empty map.
    void reserve_taxa_sorted(const std::vector<tax_t> &taxa) {
        if(!taxa_.empty() || size()) throw std::logic_error("CompactKmerMap: reserve_taxa_sorted needs an empty map");
        if(taxa.size() > max_taxa_) throw std::length_error("CompactKmerMap: more taxa than this map was sized for");
        taxa_ = taxa;
        ids_.clear();
        for(size_t i = 0; i < taxa_.size(); ++i) ids_.emplace(taxa_[i], uint32_t(i + 1));
    }
    // Drop the taxid -> dense id index (needed only by set()); it is rebuilt
    // on demand if set() is called again.
    void finalize() {decltype(ids_)().swap(ids_);}
    // f(kmer, tax) for every entry, in table order.
    template<typename F>
    void for_each(F &&f) const {
        for(auto it = t_.begin(); it != t_.end(); ++it) f(KeyT(*it), taxa_[it.mapped() - 1]);
    }
    const table_type &table() const {return t_;}
    void write(std::FILE *fp) const {
        const uint64_t nt = taxa_.size(), mt = max_taxa_;
        if(std::fwrite(&nt, sizeof(nt), 1, fp) != 1 || std::fwrite(&mt, sizeof(mt), 1, fp) != 1
           || (nt && std::fwrite(taxa_.data(), sizeof(tax_t), nt, fp) != nt))
            throw std::runtime_error("CompactKmerMap: write failed");
        t_.write(fp);
    }
    static BasicCompactKmerMap read(std::FILE *fp) {
        uint64_t nt, mt;
        if(std::fread(&nt, sizeof(nt), 1, fp) != 1 || std::fread(&mt, sizeof(mt), 1, fp) != 1)
            throw std::runtime_error("CompactKmerMap: truncated taxon table header");
        std::vector<tax_t> taxa(nt);
        if(nt && std::fread(taxa.data(), sizeof(tax_t), nt, fp) != nt)
            throw std::runtime_error("CompactKmerMap: truncated taxon table");
        BasicCompactKmerMap ret(table_type::read(fp));
        if(ret.t_.vbits() < bits_for(nt))
            throw std::runtime_error("CompactKmerMap: value width too small for the taxon table");
        ret.taxa_ = std::move(taxa);
        ret.max_taxa_ = mt;
        return ret;
    }

private:
    table_type t_;
    std::vector<tax_t> taxa_;                       // dense id - 1 -> taxid
    bns::umap<tax_t, uint32_t> ids_;       // taxid -> dense id (build-time index)
    size_t max_taxa_ = 0;
    KeyT mask_ = mask_of(t_.kbits());
    static KeyT mask_of(unsigned b) {return b >= sizeof(KeyT) * CHAR_BIT ? ~KeyT(0): ((KeyT(1) << b) - 1);}

    explicit BasicCompactKmerMap(table_type &&t): t_(std::move(t)), mask_(mask_of(t_.kbits())) {}

    // Distinct taxa of a khash_t(c) with a taxid -> dense id lookup for the
    // fill pass. Taxids below DIRECT_MAX (NCBI ids reach ~3.3M) use a direct
    // array, so the per-k-mer lookup is one array read; larger custom ids
    // fall back to a hash map. Dense ids follow sorted taxid order either way.
    struct TaxIndex {
        static constexpr tax_t DIRECT_MAX = tax_t(1) << 26;   // <= 256 MiB of u32 at the cap
        std::vector<tax_t> taxa;                             // sorted distinct taxa
        std::vector<uint32_t> direct;                        // taxid -> dense id (when max < DIRECT_MAX)
        bns::umap<tax_t, uint32_t> hashed;          // fallback
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
            bns::uset<tax_t> seen;
            for(khint_t ki = kh_begin(map); ki != kh_end(map); ++ki)
                if(kh_exist(map, ki) && seen.insert(kh_val(map, ki)).second) ix.taxa.push_back(kh_val(map, ki));
            std::sort(ix.taxa.begin(), ix.taxa.end());
            ix.hashed.reserve(ix.taxa.size());
            for(size_t i = 0; i < ix.taxa.size(); ++i) ix.hashed.emplace(ix.taxa[i], uint32_t(i + 1));
        }
        return ix;
    }
    BasicCompactKmerMap(const khash_t(c) *map, TaxIndex &&ix, double max_load, unsigned key_bits)
        : t_(q_for(kh_size(map), max_load), bits_for(ix.taxa.size()), hasher_type(key_bits)) {
        for(khint_t ki = kh_begin(map); ki != kh_end(map); ++ki) {
            if(!kh_exist(map, ki)) continue;
            if(!in_domain(kh_key(map, ki)))
                throw std::invalid_argument("CompactKmerMap: a database k-mer exceeds key_bits");
            t_.assign(kh_key(map, ki), ix.id(kh_val(map, ki)));
        }
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
using CompactKmerMap = BasicCompactKmerMap<u64>;
using CompactKmerMap128 = BasicCompactKmerMap<u128>;

// Compact classification database file: magic "BNSCKM03", then k, w and the
// spacing exactly as Database<khash_t(c)> writes them, the key encoding
// (alphabet InputType, translate flag), then the dense taxon table (count,
// capacity, taxids) and the cht_map_t stream (whose key width is the hashed
// tier's hash width). Uncompressed. "BNSCKM02" files (no key encoding field)
// read as DNA, untranslated. `bonsai compact` converts a khash database into one.
// A database stored as P compact maps, each holding the keys of one shard
// (kmer_shard). Lookups route each key to its shard; get_many keeps the order of
// its input, so classification is identical to the unsharded database. Memory is
// the sum of the shards; the point is the build (one shard at a time).
template<typename KeyT>
class ShardedKmerMap {
public:
    using map_type = BasicCompactKmerMap<KeyT>;
    using key_type = KeyT;
    static constexpr tax_t NOT_FOUND = map_type::NOT_FOUND;
    std::vector<std::unique_ptr<map_type>> shards;
    tax_t get(KeyT k) const {return shards[kmer_shard(k, shards.size())]->get(k);}
    template<typename It, typename Ot>
    void get_many(It f, It l, Ot o, std::vector<uint64_t> &scratch, size_t ahead = 16) const {
        const size_t n = std::distance(f, l), P = shards.size();
        thread_local std::vector<std::vector<KeyT>> keys;
        thread_local std::vector<std::vector<size_t>> where;
        thread_local std::vector<std::vector<tax_t>> got;
        keys.resize(P); where.resize(P); got.resize(P);
        for(size_t s = 0; s < P; ++s) keys[s].clear(), where[s].clear();
        size_t i = 0;
        for(It it = f; it != l; ++it, ++i) {
            const unsigned s = kmer_shard(KeyT(*it), P);
            keys[s].push_back(KeyT(*it));
            where[s].push_back(i);
        }
        thread_local std::vector<tax_t> out;
        out.assign(n, NOT_FOUND);
        for(size_t s = 0; s < P; ++s) {
            if(keys[s].empty()) continue;
            got[s].resize(keys[s].size());
            shards[s]->get_many(keys[s].begin(), keys[s].end(), got[s].begin(), scratch, ahead);
            for(size_t j = 0; j < keys[s].size(); ++j) out[where[s][j]] = got[s][j];
        }
        for(size_t j = 0; j < n; ++j) *o++ = out[j];
    }
    template<typename F> void for_each(F &&f) const {for(const auto &m: shards) m->for_each(f);}
    size_t size() const {size_t n = 0; for(const auto &m: shards) n += m->size(); return n;}
    size_t bytes_used() const {size_t n = 0; for(const auto &m: shards) n += m->bytes_used(); return n;}
    unsigned key_bits() const {return shards.front()->key_bits();}
    uint32_t value_bits() const {uint32_t v = 0; for(const auto &m: shards) v = std::max(v, m->value_bits()); return v;}
};

// Per-taxon k-mer counts sidecar (PI 2026-10-21): the -P/-U coverage
// denominator, computed once at database BUILD time (bonsai build) and
// serialized next to the database as <db>.taxacounts, so classify never
// rescans the whole database for it. classify loads it through
// CompactDatabase::taxon_counts() and falls back to the one O(database)
// scan only for databases built before sidecars existed (or when k or the
// canonicalization flag disagrees with the run's).
//
// File: magic "BNSKMC01" (8), u32 version, u32 k, u32 canonical, u64
// ntotal (sum of counts), u32 ntaxa, then ntaxa records of
// (u32 taxid, u64 count). Little-endian host order. Version 1.
struct TaxonCountsFile {
    static constexpr char MAGIC[9] = "BNSKMC01";
    static constexpr uint32_t VERSION = 1;
    bool ok = false;
    uint32_t k = 0;
    bool canonical = true;
    uint64_t ntotal = 0;
    bns::umap<tax_t, uint64_t> counts;
};

inline void write_taxon_counts(const char *dbpath, uint32_t k, bool canonical,
                               const bns::umap<tax_t, uint64_t> &counts) {
    const std::string path = std::string(dbpath) + ".taxacounts";
    std::FILE *fp = std::fopen(path.data(), "wb");
    if(!fp) throw std::runtime_error("could not open " + path);
    uint64_t ntotal = 0;
    for(const auto &kv: counts) ntotal += kv.second;
    const uint32_t ntaxa = uint32_t(counts.size());
    const uint32_t canon = canonical;
    if(std::fwrite(TaxonCountsFile::MAGIC, 1, 8, fp) != 8
       || std::fwrite(&TaxonCountsFile::VERSION, sizeof(uint32_t), 1, fp) != 1
       || std::fwrite(&k, sizeof(k), 1, fp) != 1
       || std::fwrite(&canon, sizeof(canon), 1, fp) != 1
       || std::fwrite(&ntotal, sizeof(ntotal), 1, fp) != 1
       || std::fwrite(&ntaxa, sizeof(ntaxa), 1, fp) != 1)
        throw std::runtime_error("could not write " + path);
    for(const auto &kv: counts) {
        const uint32_t t = kv.first;
        const uint64_t n = kv.second;
        if(std::fwrite(&t, sizeof(t), 1, fp) != 1 || std::fwrite(&n, sizeof(n), 1, fp) != 1)
            throw std::runtime_error("could not write " + path);
    }
    std::fclose(fp);
    LOG_INFO("Wrote %s: %u taxa, %llu k-mers\n", path.data(), ntaxa, (unsigned long long)ntotal);
}

inline TaxonCountsFile read_taxon_counts(const char *dbpath) {
    TaxonCountsFile ret;
    const std::string path = std::string(dbpath) + ".taxacounts";
    std::FILE *fp = std::fopen(path.data(), "rb");
    if(!fp) return ret;   // databases built before sidecars: classify falls back to the scan
    char magic[8];
    uint32_t version;
    uint32_t canon = 1, ntaxa = 0;
    try {
        if(std::fread(magic, 1, 8, fp) != 8 || std::memcmp(magic, TaxonCountsFile::MAGIC, 8)
           || std::fread(&version, sizeof(version), 1, fp) != 1 || version != TaxonCountsFile::VERSION
           || std::fread(&ret.k, sizeof(ret.k), 1, fp) != 1
           || std::fread(&canon, sizeof(canon), 1, fp) != 1
           || std::fread(&ret.ntotal, sizeof(ret.ntotal), 1, fp) != 1
           || std::fread(&ntaxa, sizeof(ntaxa), 1, fp) != 1)
            throw std::runtime_error("bad header");
        ret.canonical = canon != 0;
        for(uint32_t i = 0; i < ntaxa; ++i) {
            uint32_t t;
            uint64_t n;
            if(std::fread(&t, sizeof(t), 1, fp) != 1 || std::fread(&n, sizeof(n), 1, fp) != 1)
                throw std::runtime_error("truncated");
            ret.counts.emplace(t, n);
        }
        ret.ok = true;
    } catch(const std::exception &) {
        LOG_WARNING("%s is missing or malformed; classify rescans the database for taxon k-mer counts\n", path.data());
        ret = TaxonCountsFile{};
    }
    std::fclose(fp);
    return ret;
}

template<typename Map>
bns::umap<tax_t, uint64_t> count_taxa(const Map &m) {
    bns::umap<tax_t, uint64_t> ret;
    m.for_each([&](auto, tax_t t) {++ret[t];});
    return ret;
}

struct CompactDatabase {
    unsigned k_ = 0, w_ = 0;
    spvec_t s_;
    std::unique_ptr<CompactKmerMap> db_;        // keys of up to 64 bits
    std::unique_ptr<CompactKmerMap128> db128_;  // wider keys (DNA 32 < k <= 64 exactly)
    std::unique_ptr<ShardedKmerMap<u64>> sh_;      // sharded (BNSCKM04): keys of up to 64 bits
    std::unique_ptr<ShardedKmerMap<u128>> sh128_;  // sharded, wider keys
    InputType rht_ = DNA;                       // key alphabet
    bool translate_ = false;                    // references were 6-frame translated at build (informational)
    // Per-taxon k-mer counts (the -P/-U coverage denominator), loaded from
    // the build-time sidecar when present; see TaxonCountsFile below.
    TaxonCountsFile txc_;
    static constexpr char MAGIC[9] = "BNSCKM03";
    static constexpr char MAGIC_V2[9] = "BNSCKM02";
    static constexpr char MAGIC_V4[9] = "BNSCKM04";   // v3 header, u32 shard count, then the shard maps
    bool wide() const {return bool(db128_) || bool(sh128_);}
    bool sharded() const {return bool(sh_) || bool(sh128_);}
    // the KeySpec to classify with: reads are nucleotide, so a protein database
    // translates them (whether its references were translated or protein FASTA);
    // hash width = the stored key width, for hashed tiers
    KeySpec keyspec() const {
        KeySpec ks{rht_, rht_ != DNA, 0};
        const unsigned kb = sh128_ ? sh128_->key_bits(): sh_ ? sh_->key_bits(): wide() ? db128_->key_bits(): db_->key_bits();
        if(wide() ? KmerKeys<score::Lex, u128>::tier_for(k_, rht_) == KmerKeys<score::Lex, u128>::HASHED128
                  : KmerKeys<>::tier_for(k_, rht_) == KmerKeys<>::HASHED64) ks.hash_bits = kb;
        return ks;
    }

    static bool is_compact(const char *fn) {
        std::FILE *fp = std::fopen(fn, "rb");
        if(!fp) return false;
        char m[8];
        const bool ret = std::fread(m, 1, 8, fp) == 8 && (std::memcmp(m, MAGIC, 8) == 0 || std::memcmp(m, MAGIC_V2, 8) == 0
                                                          || std::memcmp(m, MAGIC_V4, 8) == 0);
        std::fclose(fp);
        return ret;
    }
    explicit CompactDatabase(const char *fn) {
        std::FILE *fp = std::fopen(fn, "rb");
        if(!fp) LOG_EXIT("Could not open %s for reading.\n", fn);
        char m[8];
        if(std::fread(m, 1, 8, fp) != 8 || (std::memcmp(m, MAGIC, 8) && std::memcmp(m, MAGIC_V2, 8) && std::memcmp(m, MAGIC_V4, 8)))
            throw std::runtime_error(std::string("Not a compact bonsai database: ") + fn);
        const bool v4 = std::memcmp(m, MAGIC_V4, 8) == 0;
        const bool v3 = v4 || std::memcmp(m, MAGIC, 8) == 0;
        if(std::fread(&k_, sizeof(k_), 1, fp) != 1 || std::fread(&w_, sizeof(w_), 1, fp) != 1)
            throw std::runtime_error("Error: Could not read k/w from compact database");
        s_ = spvec_t(k_ - 1);
        for(auto &v: s_) {
            uint8_t bv;
            if(std::fread(&bv, sizeof(bv), 1, fp) != 1) throw std::runtime_error("Error: Could not read spacing from file");
            v = bv;
        }
        if(v3) {
            uint32_t enc[2];
            if(std::fread(enc, sizeof(enc), 1, fp) != 1) throw std::runtime_error("Error: Could not read key encoding");
            rht_ = InputType(enc[0]);
            translate_ = enc[1] != 0;
        }
        uint32_t nshards = 1;
        if(v4 && std::fread(&nshards, sizeof(nshards), 1, fp) != 1) throw std::runtime_error("Error: Could not read shard count");
        // the key width (cht_map_t header, after the taxon table) picks the map type
        const long pos = std::ftell(fp);
        uint64_t nt, mt;
        char cm[8];
        uint32_t kb = 0;
        if(std::fread(&nt, sizeof(nt), 1, fp) != 1 || std::fread(&mt, sizeof(mt), 1, fp) != 1
           || std::fseek(fp, long(nt * sizeof(tax_t)), SEEK_CUR) != 0
           || std::fread(cm, 1, 8, fp) != 8 || std::fread(&kb, sizeof(kb), 1, fp) != 1)
            throw std::runtime_error("Error: truncated compact database");
        std::fseek(fp, pos, SEEK_SET);
        if(v4) {
            if(kb > 64) {
                sh128_.reset(new ShardedKmerMap<u128>);
                for(uint32_t i = 0; i < nshards; ++i) sh128_->shards.emplace_back(new CompactKmerMap128(CompactKmerMap128::read(fp)));
            } else {
                sh_.reset(new ShardedKmerMap<u64>);
                for(uint32_t i = 0; i < nshards; ++i) sh_->shards.emplace_back(new CompactKmerMap(CompactKmerMap::read(fp)));
            }
        } else if(kb > 64) db128_.reset(new CompactKmerMap128(CompactKmerMap128::read(fp)));
        else        db_.reset(new CompactKmerMap(CompactKmerMap::read(fp)));
        std::fclose(fp);
        txc_ = read_taxon_counts(fn);
    }
    const TaxonCountsFile &taxon_counts() const {return txc_;}
    template<typename KeyT>
    static void write(const char *fn, unsigned k, unsigned w, const spvec_t &s, const BasicCompactKmerMap<KeyT> &m,
                      const KeySpec &spec = KeySpec{}) {
        std::FILE *fp = std::fopen(fn, "wb");
        if(!fp) LOG_EXIT("Could not open %s for writing.\n", fn);
        if(std::fwrite(MAGIC, 1, 8, fp) != 8) throw std::runtime_error("Error writing database");
        if(std::fwrite(&k, sizeof(k), 1, fp) != 1 || std::fwrite(&w, sizeof(w), 1, fp) != 1)
            throw std::runtime_error("Error writing database");
        for(const auto v: s) {
            const uint8_t bv = static_cast<uint8_t>(v);
            if(std::fwrite(&bv, sizeof(bv), 1, fp) != 1) throw std::runtime_error("Error writing database");
        }
        const uint32_t enc[2] = {uint32_t(spec.rht), uint32_t(spec.translate)};
        if(std::fwrite(enc, sizeof(enc), 1, fp) != 1) throw std::runtime_error("Error writing database");
        m.write(fp);
        std::fclose(fp);
    }
    // Sharded database file: the v3 header, the shard count, then each shard's map,
    // copied from the per-shard files written during the build (so only one shard
    // is ever in memory).
    static void write_sharded(const char *fn, unsigned k, unsigned w, const spvec_t &s, const std::vector<std::string> &shard_files,
                              const KeySpec &spec = KeySpec{}) {
        std::FILE *fp = std::fopen(fn, "wb");
        if(!fp) LOG_EXIT("Could not open %s for writing.\n", fn);
        if(std::fwrite(MAGIC_V4, 1, 8, fp) != 8 || std::fwrite(&k, sizeof(k), 1, fp) != 1 || std::fwrite(&w, sizeof(w), 1, fp) != 1)
            throw std::runtime_error("Error writing database");
        for(const auto v: s) {
            const uint8_t bv = static_cast<uint8_t>(v);
            if(std::fwrite(&bv, sizeof(bv), 1, fp) != 1) throw std::runtime_error("Error writing database");
        }
        const uint32_t enc[2] = {uint32_t(spec.rht), uint32_t(spec.translate)}, n = uint32_t(shard_files.size());
        if(std::fwrite(enc, sizeof(enc), 1, fp) != 1 || std::fwrite(&n, sizeof(n), 1, fp) != 1) throw std::runtime_error("Error writing database");
        std::vector<char> buf(1 << 22);
        for(const auto &sf: shard_files) {
            std::FILE *in = std::fopen(sf.data(), "rb");
            if(!in) throw std::runtime_error("Error reading shard " + sf);
            size_t r;
            while((r = std::fread(buf.data(), 1, buf.size(), in)) > 0)
                if(std::fwrite(buf.data(), 1, r, fp) != r) throw std::runtime_error("Error writing database");
            std::fclose(in);
        }
        std::fclose(fp);
    }
};

// CompactLcaBuilder: build the k-mer -> LCA database straight into a
// CompactKmerMap, with no intermediate khash_t(c). Each genome's k-mer set is
// merged with update_lca_map's rule: absent -> the genome's taxon; present
// with a different taxon -> lca(taxon, stored).
// - Sizing: from the cardinality estimate (phase2 already computes an HLL
//   estimate), with the taxonomy's node count as the dense-id capacity.
// - Growth: when the estimate was low, the table grows to the next size by a
//   rebuild; the peak is old + new compact tables, far below khash.
// - finish(): one final rebuild at the minimum value width for the taxa
//   actually stored, with dense ids renumbered in sorted taxid order, sized
//   for the final count.
// Content (k-mer -> LCA) is independent of genome merge order; the slot
// layout is not, since linear-probing placement depends on insertion order.
template<typename KeyT = u64>
class CompactLcaBuilder {
    using map_type = BasicCompactKmerMap<KeyT>;
    std::unique_ptr<map_type> m_;
    size_t taxa_cap_;
    double max_load_;
    unsigned key_bits_;
    size_t grows_ = 0;
    unsigned shard_ = 0, nshards_ = 1;
    void grow() {
        auto next = std::make_unique<map_type>(size_t(double(m_->size() + 1) * 2), taxa_cap_, max_load_, key_bits_);
        m_->for_each([&](KeyT k, tax_t t) {next->set(k, t);});
        m_ = std::move(next);
        ++grows_;
    }
public:
    CompactLcaBuilder(size_t estimated_kmers, size_t taxonomy_nodes, unsigned key_bits = map_type::MAX_KEY_BITS,
                      double max_load = map_type::DEFAULT_MAX_LOAD):
        taxa_cap_(std::max<size_t>(taxonomy_nodes + 1, 2)), max_load_(max_load), key_bits_(key_bits) {
        m_ = std::make_unique<map_type>(size_t(double(std::max<size_t>(estimated_kmers, 16)) * 1.05),
                                        taxa_cap_, max_load_, key_bits_);
    }
    // Sharded builds: only keys of this shard are kept (kmer_shard(key, nshards) == shard).
    void set_shard(unsigned shard, unsigned nshards) {shard_ = shard; nshards_ = std::max(1u, nshards);}
    bool mine(KeyT k) const {return nshards_ == 1 || kmer_shard(k, nshards_) == shard_;}
    size_t size() const {return m_->size();}
    size_t grows() const {return grows_;}
    const map_type &map() const {return *m_;}
    // merge one k-mer of a genome with taxon taxid (update_lca_map's rule)
    void add(KeyT kmer, const khash_t(p) *tax, tax_t taxid) {
        static int warn_missing = 1;
        const tax_t cur = m_->get(kmer);
        if(cur == map_type::NOT_FOUND) {
            if(double(m_->size() + 1) > max_load_ * double(m_->table().nslots())) grow();
            m_->set(kmer, taxid);
        } else if(cur != taxid) {
            const tax_t l = lca(tax, taxid, cur);
            if(l == UINT32_C(1) && warn_missing)
                LOG_WARNING("ancestor of %u missing from taxonomy. This is not unexpected considering the issues the NCBI taxonomy has.\n", taxid), warn_missing = 0;
            m_->set(kmer, l);
        }
    }
    void update(const khash_t(all) *set, const khash_t(p) *tax, tax_t taxid) {
        for(khiter_t ki(kh_begin(set)); ki < kh_end(set); ++ki)
            if(kh_exist(set, ki) && mine(KeyT(kh_key(set, ki)))) add(KeyT(kh_key(set, ki)), tax, taxid);
    }
    void update(const std::vector<KeyT> &set, const khash_t(p) *tax, tax_t taxid) {
        for(const KeyT k: set) if(mine(k)) add(k, tax, taxid);
    }
    // Final rebuild: minimum value width, sorted dense ids, sized for size().
    std::unique_ptr<map_type> finish() {
        std::vector<tax_t> taxa;
        {
            bns::uset<tax_t> seen;
            m_->for_each([&](KeyT, tax_t t) {if(seen.insert(t).second) taxa.push_back(t);});
        }
        std::sort(taxa.begin(), taxa.end());
        auto out = std::make_unique<map_type>(m_->size(), taxa.size(), max_load_, key_bits_);
        out->reserve_taxa_sorted(taxa);
        m_->for_each([&](KeyT k, tax_t t) {out->set(k, t);});
        out->finalize();
        m_.reset();
        return out;
    }
};

// Bytes khash_t(c) holds for its buckets (keys + vals + flags), for comparison.
inline size_t khash_bytes(const khash_t(c) *map) {
    return size_t(map->n_buckets) * (sizeof(*map->keys) + sizeof(*map->vals))
         + __ac_fsize(map->n_buckets) * sizeof(khint32_t);
}

} // namespace bns

#endif /* BONSAI_KMERMAP_H__ */
