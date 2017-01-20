#ifndef TX_H__
#define TX_H__
#include <stdexcept>
#include <thread>
#include <unordered_map>
#include <functional>
#include "klib/kthread.h"

#include "lib/feature_min.h"
#include "lib/util.h"
#include "lib/counter.h"

template<typename T>
size_t get_n_occ(T *hash) {
    size_t ret(0);
    for(khiter_t ki(0); ki != kh_end(hash); ++ki) ret += !!kh_exist(hash, ki);
    return ret;
}


namespace emp {

std::string rand_string(size_t n);
uint64_t vec_popcnt(const std::string &vec);

template<typename T>
constexpr unsigned popcount(T val) {
    return __builtin_popcount(val);
}

template<>
constexpr unsigned popcount(char val) {
    return __builtin_popcount((int)val);
}

template<>
constexpr unsigned popcount(unsigned long long val) {
    return __builtin_popcountll(val);
}

template<>
constexpr unsigned popcount(unsigned long val) {
    return __builtin_popcountl(val);
}

template<typename T>
std::uint64_t vec_popcnt(T &container) {
    auto i(container.cbegin());
    std::uint64_t ret(popcount(*i));
    while(++i != container.cend()) ret += popcount(*i);
    return ret;
}

template<typename T>
uint64_t vec_popcnt(T *p, size_t l) {
    std::uint64_t ret(popcount(*p));
    while(l--) ret += popcount(*++p);
    return ret;
}

class Taxonomy {
    khash_t(p)    *tax_map_;
    khash_t(name) *name_map_;
    uint64_t       n_syn_;
    uint64_t       ceil_;
    std::string    name_;
public:
    // Textual constructor
    Taxonomy(const char *taxnodes_path, const char *name_path, const char *name="", unsigned ceil=0):
        tax_map_(build_parent_map(taxnodes_path)),
        name_map_(build_name_hash(name_path)),
        n_syn_(0),
        ceil_(ceil ? ceil: tax_map_->n_buckets << 1),
        name_(*name ? name: rand_string(20))
    {
        LOG_DEBUG("Ceil: %" PRIu64 ". n buckets: %" PRIu64 "\n", ceil_, tax_map_->n_buckets);
    }
    // Binary constructor
    Taxonomy(const char *path, unsigned ceil=0);
    ~Taxonomy() {
        kh_destroy(p, tax_map_);
        destroy_name_hash(name_map_);
    }
    void write(const char *fn) const;
    void add_node_impl(const char *node_name, const unsigned node_id, const unsigned parent);
    void add_node(const char *node_name, const unsigned parent);
    uint64_t get_syn_count() const {return n_syn_;}
    uint64_t get_syn_ceil() const {return ceil_;}
    int has_capacity() const {return n_syn_ + 1 <= ceil_;}
    bool operator==(Taxonomy &other) const;
};


struct kg_data {
    std::vector<khash_t(all) *> &core_;
    std::vector<std::string>    &paths_;
    Spacer                      &sp_;
};
void kg_helper(void *data_, long index, int tid);

class bitmap_t;

class kgset_t {
    friend bitmap_t;

    std::vector<khash_t(all) *> core_;
    std::vector<std::string>   &paths_;

public:
    void fill(std::vector<std::string> &paths, Spacer &sp, int num_threads=-1) {
        if(num_threads < 0) num_threads = std::thread::hardware_concurrency();
        kg_data data{core_, paths, sp};
        kt_for(num_threads, &kg_helper, (void *)&data, core_.size());
    }
    kgset_t(std::vector<std::string> &paths, Spacer &sp, int num_threads=-1): paths_(paths) {
        for(size_t i(0), end(paths.size()); i != end; ++i) core_.emplace_back(kh_init(all));
        fill(paths, sp, num_threads);
    }
    ~kgset_t() {
        for(auto i: core_) khash_destroy(i);
    }
    size_t size() const {return core_.size();}

    size_t weight() const {
        size_t ret(core_[0]->n_occupied);
        for(size_t i(1), e(core_.size()); i < e; ++i) ret += core_[i]->n_occupied;
        return ret;
    }
};

class bitmap_t {
    std::unordered_map<uint64_t, std::vector<std::uint64_t>> core_;
    kgset_t &set_;

    public:

    std::unordered_map<std::uint64_t, std::vector<std::uint64_t>> fill(kgset_t &set) {
        std::unordered_map<std::uint64_t, std::vector<std::uint64_t>> tmp;
        const unsigned len((set.paths_.size() + 63) >> 6);
        khash_t(all) *h;

        for(size_t i(0); i < set.core_.size(); ++i) {
            h = set.core_[i];
            for(khiter_t ki(0); ki != kh_end(h); ++ki) {
                if(kh_exist(h, ki)) {
                    auto m(tmp.find(kh_key(h, ki)));
                    if(m == tmp.end()) m = tmp.emplace(kh_key(h, ki),
                                                       std::move(std::vector<std::uint64_t>(len))).first;
                    m->second[i >> 6] |= 1u << (i & 63u);
                }
            }
        }
        return tmp;
    }

    bitmap_t(kgset_t &set): set_(set) {
        auto tmp(fill(set));
#if !NDEBUG
        size_t n_passed(0), total(tmp.size());
#endif
        unsigned bitsum;
        for(auto &i: tmp) {
            bitsum = vec_popcnt(i.second);
            if(bitsum != 1 && bitsum != set.paths_.size()) {
#if !NDEBUG
                ++n_passed;
                //LOG_DEBUG("Number passed: %zu. bitsum: %u\n", n_passed, bitsum);
#endif
                core_.emplace(i.first, i.second);
            }
        }
        LOG_DEBUG("Keeping %zu of %zu kmers for bit patterns which are not exactly compressed by the taxonomy heuristic.\n",
                  n_passed, total);
    }
    count::Counter<std::vector<std::uint64_t>> to_counter() {
        count::Counter<std::vector<std::uint64_t>> ret;
        for(auto &pair: core_) ret.add(pair.second);
        ret.set_nelem(set_.size());
        return ret;
    }
};


template<typename T>
constexpr unsigned lazy_popcnt(T val);

template<typename T>
constexpr size_t spop(T &container) {
    auto i(container.cbegin());
    std::uint64_t ret(popcount(*i));
    while(++i != container.cend()) ret += lazy_popcnt(*i);
    return ret;
}


template<typename T>
constexpr unsigned lazy_popcnt(T val) {
    unsigned ret(0);
    while(val) {
        switch(val & 0xFu) {
#ifdef SANITY
            case 1: case 2: case 4: case 8: ++ret; break;
            case 3: case 5: case 6: case 9: case 10: case 12: ret += 2; break;
            case 7: case 11: case 13: case 14: ret += 3; break;
            case 15: ret += 4; break;
#else
            case 15:                                          ++ret;
            case 7: case 11: case 13: case 14:                ++ret;
            case 3: case 5: case 6: case 9: case 10: case 12: ++ret;
            case 1: case 2: case 4: case 8:                   ++ret;
#endif
        }
        val >>= 4;
    }
    return ret;
}

static_assert(lazy_popcnt(37774) == popcount(37774), "popcnt failed");
static_assert(lazy_popcnt(3773374) == popcount(3773374), "popcnt failed");
static_assert(lazy_popcnt(0xff4cfa44) == popcount(0xff4cfa44), "popcnt failed");
static_assert(lazy_popcnt(0x1319) == popcount(0x1319), "popcnt failed");

template<typename T>
int _kh_eq(T *h1, T *h2) {
    return (h1->n_occupied == h2->n_occupied && h1->n_buckets == h2->n_buckets);
}


} // namespace emp

#endif // #ifndef TX_H__
