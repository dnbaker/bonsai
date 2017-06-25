#ifndef _TX_H__
#define _TX_H__
#include <stdexcept>
#include <thread>
#include <unordered_map>
#include <functional>
#include "klib/kthread.h"

#include "lib/feature_min.h"
#include "lib/util.h"
#include "lib/counter.h"
#include "lib/bits.h"


template<typename T>
std::size_t get_n_occ(T *hash) {
    std::size_t ret(0);
    for(khiter_t ki(0); ki != kh_end(hash); ++ki) ret += !!kh_exist(hash, ki);
    return ret;
}


namespace emp {


using namespace std::literals;

class Taxonomy {
    khash_t(p)    *tax_map_;
    khash_t(name) *name_map_;
    std::uint64_t       n_syn_;
    std::uint64_t       ceil_;
    std::string    name_;
public:
    // Textual constructor
    Taxonomy(const char *taxnodes_path, const char *name_path, const char *name="", unsigned ceil=0):
        tax_map_(build_parent_map(taxnodes_path)),
        name_map_(build_name_hash(name_path)),
        n_syn_(0),
        ceil_(ceil ? ceil: tax_map_->n_buckets << 1),
        name_(*name ? name: "Taxonomy "s + rand_string(20))
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
    std::uint64_t get_syn_count() const {return n_syn_;}
    std::uint64_t get_syn_ceil() const {return ceil_;}
    int has_capacity() const {return n_syn_ + 1 <= ceil_;}
    bool operator==(Taxonomy &other) const;
    operator khash_t(p) *() {
        return tax_map_;
    }
    operator khash_t(name) *() {
        return name_map_;
    }
};


struct kg_data {
    std::vector<khash_t(all) *> &core_;
    std::vector<std::string>    &paths_;
    const Spacer                &sp_;
    const khash_t(all)          *acceptable_;
};

struct kg_list_data {
    std::vector<khash_t(all) *> &core_;
    std::vector<std::forward_list<std::string>*> &fl_;
    const Spacer                &sp_;
    const khash_t(all)          *acceptable_;
};

void kg_helper(void *data_, long index, int tid);
void kg_list_helper(void *data_, long index, int tid);

class kgset_t {

    std::vector<khash_t(all) *> core_;
    std::vector<std::string>    paths_;
    const khash_t(all)         *acceptable_;
    std::unordered_map<std::uint32_t, std::forward_list<std::string>> *fl_;
    std::vector<tax_t>         taxes_;

public:
    void fill(std::vector<std::string> &paths, const Spacer &sp, int num_threads=-1) {
        if(num_threads < 0) num_threads = std::thread::hardware_concurrency();
        kg_data data{core_, paths, sp, acceptable_};
        kt_for(num_threads, &kg_helper, (void *)&data, core_.size());
    }
    void fill(std::unordered_map<std::uint32_t, std::forward_list<std::string>> *path_map, const Spacer &sp, int num_threads=-1) {
        auto &pm(*path_map);
        if(num_threads < 0) num_threads = std::thread::hardware_concurrency();
        std::vector<std::forward_list<std::string>*> tmpfl;
        taxes_.reserve(pm.size());
        LOG_DEBUG("Tax reserved size %zu\n", taxes_.capacity());
        for(auto &pair: pm) {
            taxes_.push_back(pair.first);
            assert(emp::size(pair.second) > 0);
            tmpfl.push_back(&pair.second);
        }
        LOG_DEBUG("Tax filled size %zu. Core size: %zu\n", taxes_.size(), core_.size());
        kg_list_data data{core_, tmpfl, sp, acceptable_};
        kt_for(num_threads, &kg_list_helper, (void *)&data, core_.size());
    }
    auto &get_taxes() {return taxes_;} // Who'd want that?
    const std::vector<khash_t(all) *> &get_core() const {return core_;}
    std::vector<std::string>    &get_paths() {return paths_;}
    // Encoding constructor
    kgset_t(typename std::vector<std::string>::iterator begin, typename std::vector<std::string>::iterator end,
            const Spacer &sp, int num_threads=-1, const khash_t(all) *acc=nullptr): paths_(begin, end), acceptable_(acc), fl_(nullptr) {
        core_.reserve(end - begin);
        for(std::size_t i(0), end(paths_.size()); i != end; ++i) core_.emplace_back(kh_init(all));
        fill(paths_, sp, num_threads);
    }
    // Encoding constructor
    kgset_t(std::vector<std::string> &paths, const Spacer &sp, int num_threads=-1, const khash_t(all) *acc=nullptr):
        kgset_t(std::begin(paths), std::end(paths), sp, num_threads, acc) {
    }
    kgset_t(std::unordered_map<std::uint32_t, std::forward_list<std::string>> &list,
            const Spacer &sp, int num_threads=-1, const khash_t(all) *acc=nullptr): acceptable_(acc), fl_(&list) {
        LOG_DEBUG("Acc? %p\n", (void *)acc);
        if(list.size() == 0) LOG_EXIT("List size is 0\n");
        core_.reserve(list.size());
        while(core_.size() < list.size()) core_.emplace_back(kh_init(all));
        assert(core_.size() > 0);
        fill(fl_, sp, num_threads);
    }
    ~kgset_t() {
        for(auto i: core_) khash_destroy(i);
    }
    std::size_t size() const {return core_.size();}

    std::size_t weight() const {
        if(!size()) return 0;
        std::size_t ret(kh_size(core_[0]));
        for(std::size_t i(1), e(core_.size()); i < e; ++i) ret += kh_size(core_[i]);
        return ret;
    }
};


template<typename T>
constexpr std::size_t spop(T &container) {
    assert(container.size());
    auto i(container.cbegin());
    std::uint64_t ret(popcnt::popcount(*i));
    while(++i != container.cend()) ret += popcnt::popcount(*i);
    return ret;
}

#if 0

template<typename T>
constexpr unsigned lazy_popcnt(T val) {
    unsigned ret(0);
    while(val) {
        switch(val & 0xFu) {
#ifdef SANITY
            case 1: case 2: case 4: case 8:                   ++ret;    break;
            case 3: case 5: case 6: case 9: case 10: case 12: ret += 2; break;
            case 7: case 11: case 13: case 14:                ret += 3; break;
            case 15:                                          ret += 4; break;
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

static_assert(lazy_popcnt(37774) == popcnt::popcount(37774), "popcnt failed");
static_assert(lazy_popcnt(3773374) == popcnt::popcount(3773374), "popcnt failed");
static_assert(lazy_popcnt(0xff4cfa44) == popcnt::popcount(0xff4cfa44), "popcnt failed");
static_assert(lazy_popcnt(0x1319) == popcnt::popcount(0x1319), "popcnt failed");
#endif

template<typename T>
int _kh_eq(T *h1, T *h2) {
    return (h1->n_occupied == h2->n_occupied && h1->n_buckets == h2->n_buckets);
}


} // namespace emp

#endif // #ifndef _TX_H__
