#pragma once
#include "feature_min.h"
#include "counter.h"

namespace bns {

struct kg_data {
    std::vector<khash_t(all)>   &core_;
    std::vector<std::string>   &paths_;
    const Spacer                  &sp_;
    const khash_t(all)    *acceptable_;
    const bool                  canon_;
};

struct kg_list_data {
    std::vector<khash_t(all)>   &core_;
    const std::vector<const std::forward_list<std::string>*> &fl_;
    const Spacer                  &sp_;
    const khash_t(all)    *acceptable_;
    const bool                  canon_;
};

static void kg_helper(void *data_, long index, int) {
    kg_data *data(reinterpret_cast<kg_data *>(data_));
    khash_t(all) *hash(&data->core_[index]);
    int khr;
    Encoder<score::Lex> enc(data->sp_, data->canon_);
    enc.for_each([&](u64 min) {
        //LOG_INFO("Kmer is %s\n", data->sp_.to_string(min).data());
        if(!data->acceptable_ || (kh_get(all, data->acceptable_, min) != kh_end(data->acceptable_)))
            kh_put(all, hash, min, &khr);
        assert(kh_size(hash));
    }, data->paths_[index].data());
}

static void kg_list_helper(void *data_, long index, int) {
    kg_list_data &data(*reinterpret_cast<kg_list_data *>(data_));
    auto &list(*data.fl_[index]);
    khash_t(all) *hash(&data.core_[index]);
    int khr;
    Encoder<score::Lex> enc(data.sp_, data.canon_);
    enc.for_each([&](u64 min) {
        if(!data.acceptable_ || (kh_get(all, data.acceptable_, min) != kh_end(data.acceptable_)))
            khash_put(hash, min, &khr);
    }, list);
}

class kgset_t {

    std::vector<khash_t(all)>    core_;
    std::vector<std::string>    paths_;
    const khash_t(all)         *acceptable_;
    const std::unordered_map<u32, std::forward_list<std::string>> *fl_;
    std::vector<tax_t>          taxes_;

public:
    void fill(std::vector<std::string> &paths, const Spacer &sp, bool canonicalize=true, int num_threads=-1) {
        if(num_threads < 0) num_threads = std::thread::hardware_concurrency();
        kg_data data{core_, paths, sp, acceptable_, canonicalize};
        ForPool pool(num_threads);
        pool.forpool(&kg_helper, reinterpret_cast<void *>(&data), core_.size());
    }
    void fill(const std::unordered_map<u32, std::forward_list<std::string>> *path_map, const Spacer &sp, bool canonicalize=true, int num_threads=-1) {
        auto &pm(*path_map);
        if(num_threads < 0) num_threads = std::thread::hardware_concurrency();
        std::vector<const std::forward_list<std::string>*> tmpfl;
        taxes_.clear(), taxes_.reserve(std::max(pm.size(), size_t(4)));
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
            assert(::bns::size(pair.second) > 0);
            tmpfl.push_back(&pair.second);
        }
        kg_list_data data{core_, tmpfl, sp, acceptable_, canonicalize};
        ForPool pool(num_threads);
        pool.forpool(&kg_list_helper, reinterpret_cast<void *>(&data), core_.size());
    }
    const auto &get_taxes()                        const {return taxes_;} // Who'd want that?
    const std::vector<khash_t(all)> &core()        const {return core_;}
    const std::vector<std::string> &paths()        const {return paths_;}
    // Encoding constructors
    kgset_t(typename std::vector<std::string>::const_iterator begin, typename std::vector<std::string>::const_iterator end,
            const Spacer &sp, bool canonicalize=true, int num_threads=-1, const khash_t(all) *acc=nullptr): paths_(begin, end), acceptable_(acc), fl_(nullptr) {
        core_.reserve(end - begin);
        for(size_t i(0), end(paths_.size()); i != end; ++i) core_.emplace_back(khash_t(all){0,0,0,0,0,0,0});
        fill(paths_, sp, canonicalize, num_threads);
    }
    kgset_t(const std::vector<std::string> &paths, const Spacer &sp, bool canonicalize=true, int num_threads=-1, const khash_t(all) *acc=nullptr):
        kgset_t(std::begin(paths),
                std::end(paths), sp, canonicalize, num_threads, acc) {}

    // Forward list constructor
    kgset_t(const std::unordered_map<u32, std::forward_list<std::string>> &list,
            const Spacer &sp, bool canonicalize=true, int num_threads=-1, const khash_t(all) *acc=nullptr): acceptable_(acc), fl_(&list) {
        if(list.size() == 0) LOG_EXIT("List size is 0\n");
        core_.reserve(list.size());
        while(core_.size() < list.size()) core_.emplace_back(khash_t(all){0,0,0,0,0,0,0});
        assert(core_.size());
        fill(fl_, sp, canonicalize, num_threads);
    }

    ~kgset_t() {
        for(auto &i: core_) std::free(i.keys), std::free(i.vals), std::free(i.flags);
    }
    size_t size() const {return core_.size();}

    size_t weight() const {
        if(!size()) return 0;
        size_t ret(kh_size(&core_[0]));
        for(size_t i(1), e(core_.size()); i < e; ++i) ret += kh_size(&core_[i]);
        return ret;
    }
    void print_weights(std::FILE *fp=stderr) const {
        //for(const auto &kh: core_)
            //std::fprintf(fp, "Occupancy of %zu\n", kh_size(&kh));
    }
};

} // namespace bns
