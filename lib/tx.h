#ifndef TX_H__
#define TX_H__
#include <stdexcept>
#include <thread>
#include <unordered_map>
#include "klib/kthread.h"

#include "lib/feature_min.h"
#include "lib/util.h"
#include "lib/counter.h"

namespace emp {

template<typename T>
int _kh_eq(T *h1, T *h2) {
    print_khash(h1);
    print_khash(h2);
    LOG_DEBUG("Is %s\n", (h1->n_occupied == h2->n_occupied && h1->n_buckets == h2->n_buckets) ?
                         "true": "false");
    return (h1->n_occupied == h2->n_occupied && h1->n_buckets == h2->n_buckets);
}

std::string rand_string(size_t n);

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

template<typename T>
unsigned popcount(T val);
template<>
unsigned popcount(unsigned long val);
template<>
unsigned popcount(unsigned long long val);
uint64_t vec_popcnt(std::string &vec);

class bitmap_t {
    std::unordered_map<uint64_t, std::string> core_;
    kgset_t &set_;

    std::unordered_map<uint64_t, std::string> fill(kgset_t &set) {
        std::unordered_map<uint64_t, std::string> tmp;
        const unsigned len((set.paths_.size() + CHAR_BIT - 1) / CHAR_BIT);
        khash_t(all) *h;

        for(size_t i(0); i < set.core_.size(); ++i) {
            h = set.core_[i];
            for(khiter_t ki(0); ki != kh_end(h); ++ki) {
                if(kh_exist(h, ki)) {
                    auto m(tmp.find(kh_key(h, ki)));
                    if(m == tmp.end()) m = tmp.emplace(kh_key(h, ki),
                                                       std::move(std::string(len, '\0'))).first;
                    m->second[i >> 3] |= 1 << (i & 7);
                }
            }
        }
        return tmp;
    }

    bitmap_t(kgset_t &set): set_(set) {
        auto tmp(fill(set));
        for(auto &i: tmp) {
            const unsigned bitsum(vec_popcnt(i.second));
            if(bitsum == 1 || bitsum == set.paths_.size()) continue;
            core_.emplace(i.first, i.second);
        }
    }
    count::Counter<std::string> to_counter() {
        count::Counter<std::string> ret;
        for(auto &pair: core_) ret.add(pair.second);
        return ret;
    }
};

} // namespace emp

#endif // #ifndef TX_H__
