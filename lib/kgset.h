#pragma once
#include "lib/feature_min.h"
#include "lib/counter.h"

namespace emp {

struct kg_data {
    std::vector<khash_t(all) *> &core_;
    std::vector<std::string>    &paths_;
    const Spacer                &sp_;
    const khash_t(all)          *acceptable_;
};

struct kg_list_data {
    std::vector<khash_t(all) *> &core_;
    const std::vector<const std::forward_list<std::string>*> &fl_;
    const Spacer                &sp_;
    const khash_t(all)          *acceptable_;
};

void kg_helper(void *data_, long index, int tid);
void kg_list_helper(void *data_, long index, int tid);
class kgset_t {

    std::vector<khash_t(all) *> core_;
    std::vector<std::string>    paths_;
    const khash_t(all)         *acceptable_;
    const std::unordered_map<u32, std::forward_list<std::string>> *fl_;
    std::vector<tax_t>         taxes_;

public:
    void fill(std::vector<std::string> &paths, const Spacer &sp, int num_threads=-1) {
        if(num_threads < 0) num_threads = std::thread::hardware_concurrency();
        kg_data data{core_, paths, sp, acceptable_};
        kt_for(num_threads, &kg_helper, (void *)&data, core_.size());
    }
    void fill(const std::unordered_map<u32, std::forward_list<std::string>> *path_map, const Spacer &sp, int num_threads=-1) {
        auto &pm(*path_map);
        if(num_threads < 0) num_threads = std::thread::hardware_concurrency();
        std::vector<const std::forward_list<std::string>*> tmpfl;
        taxes_.clear(), taxes_.reserve(std::max(pm.size(), (size_t)4));
        LOG_DEBUG("Tax reserved size %zu\n", taxes_.capacity());
#if !NDEBUG
        for(const auto &pair: pm) {
            std::cerr << "Tax: " << pair.first;
            std::cerr << "Paths: [";
            for(const auto &el: pair.second) std::cerr << el << ", ";
            std::cerr << "].\n";
        }
#endif
        for(const auto &pair: pm) {
            taxes_.push_back(pair.first);
            assert(emp::size(pair.second) > 0);
            tmpfl.push_back(&pair.second);
        }
        LOG_DEBUG("Tax filled size %zu. Core size: %zu\n", taxes_.size(), core_.size());
        kg_list_data data{core_, tmpfl, sp, acceptable_};
        kt_for(num_threads, &kg_list_helper, (void *)&data, core_.size());
    }
    const auto &get_taxes()                        const {return taxes_;} // Who'd want that?
    const std::vector<khash_t(all) *> &get_core()  const {return core_;}
    const std::vector<std::string>    &get_paths() const {return paths_;}
    // Encoding constructors
    kgset_t(typename std::vector<std::string>::const_iterator begin, typename std::vector<std::string>::const_iterator end,
            const Spacer &sp, int num_threads=-1, const khash_t(all) *acc=nullptr): paths_(begin, end), acceptable_(acc), fl_(nullptr) {
        core_.reserve(end - begin);
        for(size_t i(0), end(paths_.size()); i != end; ++i) core_.emplace_back(kh_init(all));
        fill(paths_, sp, num_threads);
    }
    kgset_t(const std::vector<std::string> &paths, const Spacer &sp, int num_threads=-1, const khash_t(all) *acc=nullptr):
        kgset_t(std::begin(paths),
                std::end(paths), sp, num_threads, acc) {}

    // Forward list constructor
    kgset_t(const std::unordered_map<u32, std::forward_list<std::string>> &list,
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
    size_t size() const {return core_.size();}

    size_t weight() const {
        if(!size()) return 0;
        size_t ret(kh_size(core_[0]));
        for(size_t i(1), e(core_.size()); i < e; ++i) ret += kh_size(core_[i]);
        return ret;
    }
};


} // namespace emp
