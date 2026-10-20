#ifndef BONSAI_TAXMAP_H__
#define BONSAI_TAXMAP_H__
// CompactTaxMap: the taxonomic tree (tax -> parent) in sketch's compact hash
// table (hll/include/sketch/cht.h) instead of khash_t(p).
//
// ENCODING: the tree is stored as parent POINTERS IN SLOT SPACE. Each taxid
// lives in one cht slot, and the slot's value is (the parent's slot + 1).
// - Keys: the 32-bit taxid is mixed by a bijection (murmur3 fmix32). The top q
//   bits pick the home slot and only the remaining 32 - q bits are stored.
// - Values: VAL_BITS = q_max + 2 bits (26) instead of a 32-bit taxid.
// - Walking up the tree (lca, node_depth, node_dist) costs one hash per QUERY
//   taxid. Every further step is a packed-array read, not a hash-and-probe.
// - Per slot: 32 - q key bits + 26 value bits + 1 occupancy bit + 8
//   displacement bits (used for exact matching and key decoding). That is
//   about 45 bits at NCBI scale (q = 22), against khash_t(p)'s 66 bits per
//   bucket at a lower load.
// - Slots are durable handles: linear probing without deletion never moves
//   an entry.
//
// SPECIAL VALUES:
//   ROOT_PARENT: the parent is 0 (the root's parent, or any row with parent 0).
//   DANGLING: the parent taxid is not (yet) in the map. It is kept in a side
//     table and re-linked automatically when that taxid is inserted, so the
//     result never depends on insertion order: parent() and the walks behave
//     exactly as khash_t(p) does for the same rows. After a consistent build
//     the side table is empty.
// A parent of tax_t(-1) is not storable (set_parent throws); build_parent_map
// already drops the malformed rows that produced it.
//
// EXACTNESS: fmix32 is a bijection, the full remainder is stored, and every
// lookup uses cht's home-exact find_slot (remainder AND displacement match).
// So there are no merged keys or false positives, as with khash. (cht's
// count() alone is not exact at a 32-bit key width.)
//
// The free functions (get_parent, node_depth, node_dist, lca) mirror the
// khash_t(p) overloads in util.h, including their missing-taxid messages and
// return values, so call sites can migrate by changing the map type.
// Regression: test/regress_compact_taxmap.cpp. Memory/speed: bench/taxmap_bench.cpp.
#include "util.h"
#include "sketch/cht.h"

#include <cmath>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace bns {

namespace taxmap_detail {
constexpr uint32_t inv_mul(uint32_t a) {  // a odd: Newton iteration for a^-1 mod 2^32
    uint32_t x = a;
    for(int i = 0; i < 5; ++i) x *= 2u - a * x;
    return x;
}
} // namespace taxmap_detail

// murmur3 fmix32 and its inverse: a bijection on [0, 2^32) that spreads the
// dense, clustered small integers taxids are.
struct TaxIdRevHash {
    static constexpr size_t key_bits = 32;
    static constexpr uint32_t C1 = 0x85ebca6bu, C2 = 0xc2b2ae35u;
    static constexpr uint32_t C1I = taxmap_detail::inv_mul(C1), C2I = taxmap_detail::inv_mul(C2);
    static_assert(uint32_t(C1 * C1I) == 1u && uint32_t(C2 * C2I) == 1u, "modular inverses");
    uint64_t hash(uint64_t key) const {
        uint32_t x = uint32_t(key);
        x ^= x >> 16; x *= C1; x ^= x >> 13; x *= C2; x ^= x >> 16;
        return x;
    }
    uint64_t invert(uint64_t h) const {
        uint32_t x = uint32_t(h);
        x ^= x >> 16;
        x *= C2I;
        x ^= (x >> 13) ^ (x >> 26);   // inverse of x ^= x >> 13 on 32 bits
        x *= C1I;
        x ^= x >> 16;
        return x;
    }
};

class CompactTaxMap {
public:
    static constexpr size_t VAL_BITS = 26;
    using table_type = sketch::cht::cht_t<uint64_t, TaxIdRevHash, VAL_BITS>;
    static constexpr tax_t NOT_FOUND = std::numeric_limits<tax_t>::max();
    static constexpr double DEFAULT_MAX_LOAD = 0.85;
    static constexpr uint64_t ROOT_PARENT = table_type::SAT;       // parent is 0
    static constexpr uint64_t DANGLING = table_type::SAT - 1;      // parent taxid not in the map
    static constexpr uint32_t MAX_Q = 24;                          // slot + 1 < DANGLING
    static_assert((uint64_t(1) << MAX_Q) < DANGLING, "slot pointers must stay below the special values");

    // Smallest table that holds n entries at load <= max_load.
    static uint32_t q_for(size_t n, double max_load = DEFAULT_MAX_LOAD) {
        uint32_t q = 2;
        while(q < MAX_Q && double(n) > max_load * std::ldexp(1., int(q))) ++q;
        if(double(n) > max_load * std::ldexp(1., int(q)))
            throw std::length_error("CompactTaxMap: more taxa than a 2^24-slot table holds at this load");
        return q;
    }

    explicit CompactTaxMap(size_t expected_size, double max_load = DEFAULT_MAX_LOAD)
        : t_(q_for(expected_size, max_load)) {}
    // Conversion from an existing khash_t(p) (same entries, same parents).
    explicit CompactTaxMap(const khash_t(p) *map, double max_load = DEFAULT_MAX_LOAD)
        : t_(q_for(kh_size(map), max_load)) {
        for(khint_t ki = kh_begin(map); ki != kh_end(map); ++ki)
            if(kh_exist(map, ki)) set_parent(kh_key(map, ki), kh_val(map, ki));
    }

    // Insert child, or overwrite its parent (kh_put + kh_val semantics, any order).
    void set_parent(tax_t child, tax_t parent) {
        if(parent == NOT_FOUND)
            throw std::invalid_argument("CompactTaxMap: parent tax_t(-1) is not storable");
        uint64_t cs = t_.find_slot(child);
        if(cs == t_.nslots()) {                       // new taxid: insert, then adopt waiting children
            cs = t_.assign(child, ROOT_PARENT);       // placeholder; overwritten below
            auto w = waiting_.find(child);
            if(w != waiting_.end()) {
                for(const uint64_t c: w->second) {
                    t_.set_value_at(c, cs + 1);
                    dangling_.erase(c);
                }
                waiting_.erase(w);
            }
        } else {
            unlink_dangling(cs);
        }
        t_.set_value_at(cs, encode_parent(cs, parent));
    }
    // Parent of `child`, or NOT_FOUND (tax_t(-1)) when absent, like get_parent(khash).
    tax_t parent(tax_t child) const {
        const uint64_t s = t_.find_slot(child);
        return s == t_.nslots() ? NOT_FOUND: parent_of_slot(s);
    }
    bool contains(tax_t child) const {return t_.find_slot(child) != t_.nslots();}
    size_t size() const {return t_.size();}
    double load_factor() const {return t_.load_factor();}
    // Table bytes plus the (normally empty) dangling side tables, approximately.
    size_t bytes_used() const {
        size_t w = 0;
        for(const auto &kv: waiting_) w += sizeof(kv) + kv.second.capacity() * sizeof(uint64_t);
        return t_.bytes_used() + dangling_.size() * (sizeof(uint64_t) + sizeof(tax_t) + 2 * sizeof(void *)) + w;
    }
    size_t num_dangling() const {return dangling_.size();}
    // f(child, parent) for every entry, in table order.
    template<typename F>
    void for_each(F &&f) const {
        for(auto it = t_.begin(); it != t_.end(); ++it) f(tax_t(*it), parent_of_slot(it.slot()));
    }
    const table_type &table() const {return t_;}

    // ---- slot-space walking (used by the free functions below) ----
    uint64_t nslots() const {return t_.nslots();}
    uint64_t slot_of(tax_t t) const {return t_.find_slot(t);}
    tax_t taxid_at(uint64_t s) const {return tax_t(t_.decode(s));}
    // Next slot up the tree from slot s. Returns nslots() when s's parent is 0
    // (top reached). When the parent taxid is absent from the map, returns
    // nslots() + 1 and sets *missing to that taxid.
    uint64_t parent_slot(uint64_t s, tax_t *missing) const {
        const uint64_t v = t_.value_at(s);
        if(v == ROOT_PARENT) return t_.nslots();
        if(v == DANGLING) {*missing = dangling_.at(s); return t_.nslots() + 1;}
        return v - 1;
    }

private:
    table_type t_;
    std::unordered_map<uint64_t, tax_t> dangling_;                 // child slot -> absent parent taxid
    std::unordered_map<tax_t, std::vector<uint64_t>> waiting_;     // absent parent taxid -> child slots

    uint64_t encode_parent(uint64_t cs, tax_t parent) {
        if(parent == 0) return ROOT_PARENT;
        const uint64_t ps = t_.find_slot(parent);
        if(ps != t_.nslots()) return ps + 1;
        dangling_[cs] = parent;
        waiting_[parent].push_back(cs);
        return DANGLING;
    }
    void unlink_dangling(uint64_t cs) {
        auto d = dangling_.find(cs);
        if(d == dangling_.end()) return;
        auto w = waiting_.find(d->second);
        if(w != waiting_.end()) {
            auto &v = w->second;
            v.erase(std::remove(v.begin(), v.end(), cs), v.end());
            if(v.empty()) waiting_.erase(w);
        }
        dangling_.erase(d);
    }
    tax_t parent_of_slot(uint64_t s) const {
        const uint64_t v = t_.value_at(s);
        if(v == ROOT_PARENT) return 0;
        if(v == DANGLING) return dangling_.at(s);
        return tax_t(t_.decode(v - 1));
    }
};

// Bytes khash_t(p) holds for its buckets (keys + vals + flags), for comparison.
inline tax_t tax_parent(const CompactTaxMap &m, tax_t t) noexcept {return m.parent(t);}

inline size_t khash_bytes(const khash_t(p) *map) {
    return size_t(map->n_buckets) * (sizeof(*map->keys) + sizeof(*map->vals))
         + __ac_fsize(map->n_buckets) * sizeof(khint32_t);
}

INLINE tax_t get_parent(const CompactTaxMap &taxmap, tax_t key) noexcept {
    return taxmap.parent(key);
}

static unsigned node_dist(const CompactTaxMap &map, tax_t leaf, tax_t root) noexcept {
    unsigned ret(0);
    if(!leaf) LOG_EXIT("leaf %u is not a child of root %u\n", leaf, root);
    uint64_t s = map.slot_of(leaf);
    const uint64_t top = map.nslots();
    for(;;) {
        if(s >= top) {   // the current taxid (leaf) is absent
            std::fprintf(stderr, "Tax ID %u missing. Abort!\n", leaf);
            std::exit(1);
        }
        ++ret;
        tax_t missing = 0;
        const uint64_t ps = map.parent_slot(s, &missing);
        leaf = ps == top ? 0: ps == top + 1 ? missing: map.taxid_at(ps);
        if(leaf == root) return ret;
        if(!leaf) break;
        s = ps == top + 1 ? top: ps;
    }
    LOG_EXIT("leaf %u is not a child of root %u\n", leaf, root);
    return ret;
}

static unsigned node_depth(const CompactTaxMap &map, tax_t a) noexcept {
    unsigned ret(0);
    if(!a) return 0;
    uint64_t s = map.slot_of(a);
    const uint64_t top = map.nslots();
    for(;;) {
        if(s >= top) {
            std::fprintf(stderr, "Tax ID %u missing. Abort!\n", a);
            std::exit(1);
        }
        ++ret;
        tax_t missing = 0;
        const uint64_t ps = map.parent_slot(s, &missing);
        if(ps == top) return ret;
        if(ps == top + 1) {a = missing; s = top; continue;}
        s = ps;
    }
}

// Same algorithm and return values as lca(const khash_t(p) *, a, b), walking
// slot pointers: the ancestors of a are collected as slots, then b climbs
// until it meets one.
static tax_t lca(const CompactTaxMap &map, tax_t a, tax_t b) noexcept {
    if(a == b) return a;
    if(b == 0) return a;
    if(a == 0) return b;
    const uint64_t top = map.nslots();
    linear::set<uint64_t> nodes;
    uint64_t s = map.slot_of(a);
    for(;;) {
        if(s >= top) {
            std::fprintf(stderr, "Missing taxid %u. Returning node b (%u)!\n", a, b);
            return (tax_t)-1;
        }
        nodes.insert(s);
        tax_t missing = 0;
        const uint64_t ps = map.parent_slot(s, &missing);
        if(ps == top) break;
        if(ps == top + 1) {a = missing; s = top; continue;}
        s = ps;
    }
    a = 0;   // as in the khash version, whose first loop ends with a == 0 (read by the message below)
    s = map.slot_of(b);
    for(;;) {
        if(s >= top) {
            std::fprintf(stderr, "Missing taxid %u. Returning node b (%u)!\n", b, a);
            return (tax_t)-1;
        }
        if(nodes.find(s) != nodes.end()) return map.taxid_at(s);
        tax_t missing = 0;
        const uint64_t ps = map.parent_slot(s, &missing);
        if(ps == top) break;
        if(ps == top + 1) {b = missing; s = top; continue;}   // b is only read by the message
        s = ps;
    }
    return 1;
}

template<typename Container, typename=typename std::enable_if<std::is_same<typename Container::value_type, tax_t>::value>::type>
tax_t lca(const CompactTaxMap &taxmap, const Container &v) noexcept {
    if(v.size() == 0) {
        fprintf(stderr, "Warning: no elements provided. Returning 0 for lca.\n");
        return 0;
    }
    auto it(v.begin());
    auto ret(*it);
    while(++it != v.end()) ret = lca(taxmap, ret, *it);
    return ret;
}

// build_parent_map's parsing rules exactly (comment/blank rows skipped,
// pipeless rows warned and skipped, later rows overwrite earlier ones, the
// root 1 forced to parent 0), but straight into a CompactTaxMap: one pass to
// size the table, one to fill it, and no intermediate khash.
inline CompactTaxMap build_compact_parent_map(const char *fn, double max_load = CompactTaxMap::DEFAULT_MAX_LOAD) {
    size_t nrows = 1;   // + the forced root
    {
        std::ifstream is(fn);
        if(!is) RUNTIME_ERROR(std::string("Failed to open taxonomy ") + fn);
        std::string line;
        while(std::getline(is, line))
            if(!line.empty() && line[0] != '#' && line[0] != '\n' && std::strchr(line.data(), '|')) ++nrows;
    }
    CompactTaxMap ret(nrows, max_load);
    std::ifstream is(fn);
    std::string line;
    while(std::getline(is, line)) {
        switch(line[0]) case '\n': case '\0': case '#': continue;
        const char *p = std::strchr(line.data(), '|');
        if(p == nullptr) {
            LOG_WARNING("Malformed line (no pipe): %s", line.data());
            continue;
        }
        ret.set_parent(std::atoi(line.data()), std::atoi(p + 1));
    }
    ret.set_parent(1, 0); // Root of the tree.
    if(ret.size() < 2) RUNTIME_ERROR(std::string("Failed to create taxmap from ") + fn);
    LOG_DEBUG("Built compact parent map of size %zu (%zu bytes) from path %s\n", ret.size(), ret.bytes_used(), fn);
    return ret;
}

} // namespace bns

#endif /* BONSAI_TAXMAP_H__ */
