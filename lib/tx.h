#ifndef TX_H__
#define TX_H__
#include <stdexcept>
#include <thread>
#include <unordered_map>
#include "klib/kthread.h"

#include "lib/feature_min.h"
#include "lib/util.h"

namespace emp {


class Taxonomy {
    khash_t(p)    *tax_map_;
    khash_t(name) *name_map_;
    uint64_t       n_syn_;
    uint64_t       ceil_;
public:
    // Textual constructor
    Taxonomy(const char *taxnodes_path, const char *name_path, unsigned ceil=0):
        tax_map_(build_parent_map(taxnodes_path)),
        name_map_(build_name_hash(name_path)),
        n_syn_(0),
        ceil_(ceil ? ceil: tax_map_->n_buckets << 1)
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
    int has_capacity() {return n_syn_ + 1 <= ceil_;}
    bool operator==(Taxonomy &other) {
        if(n_syn_ != other.n_syn_) return false;
        if(ceil_ != other.ceil_) return false;
        for(khiter_t ki(0); ki != kh_end(tax_map_); ++ki) {
            // Make sure that all keys are shared between maps and all values are equal.
        }
    }
};


namespace {
struct kg_data {
    std::vector<khash_t(all) *> &core_;
    std::vector<std::string>    &paths_;
    Spacer                      &sp_;
};
void kg_helper(void *data_, long index, int tid);
}

class bitmap_t;

class kgset_t {
    friend bitmap_t;

    std::vector<khash_t(all) *> core_;
    std::vector<std::string>   &paths_;

public:
    kgset_t(std::vector<std::string> &paths): paths_(paths) {
        for(size_t i(0), end(paths.size()); i != end; ++i) core_.emplace_back(kh_init(all));
    }
    void fill(std::vector<std::string> &paths, Spacer &sp, int num_threads=-1) {
        if(num_threads < 0) num_threads = std::thread::hardware_concurrency();
        kg_data data{core_, paths, sp};
        kt_for(num_threads, &kg_helper, (void *)&data, core_.size());
    }
    ~kgset_t() {
        for(auto i: core_) khash_destroy(i);
    }
};

template<typename T>
unsigned popcount(T val);
template<>
unsigned popcount(unsigned long val);
template<>
unsigned popcount(unsigned long long val);
uint64_t vec_popcnt(std::vector<uint64_t> &vec);

class bitmap_t {
    std::unordered_map<uint64_t, std::vector<uint64_t>> core_;
    kgset_t &set_;

    std::unordered_map<uint64_t, std::vector<uint64_t>> fill(kgset_t &set) {
        std::unordered_map<uint64_t, std::vector<uint64_t>> tmp;
        const unsigned len((set.paths_.size() + 63) / 64);
        khash_t(all) *h;

        for(size_t i(0); i < set.core_.size(); ++i) {
            h = set.core_[i];
            for(khiter_t ki(0); ki != kh_end(h); ++ki) {
                if(kh_exist(h, ki)) {
                    auto m(tmp.find(kh_key(h, ki)));
                    if(m == tmp.end())
                        m = tmp.emplace(kh_key(h, ki), std::move(std::vector<uint64_t>(len))).first;
                    m->second[i >> 6] |= 1 << (i & 63);
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
};

} // namespace emp

#endif // #ifndef TX_H__
