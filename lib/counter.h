#ifndef COUNTER_H__
#define COUNTER_H__
#include <algorithm>
#include <typeinfo>
#include <cstddef>
#include <cstdlib>
#include <climits>
#include <memory>
#include <set>
#include <unordered_map>
#include <vector>
#include "lib/logutil.h"
#include "lib/ks.h"
#include "clhash/include/clhash.h"

class rand_holder {
    void *random_;
public:
    rand_holder(): random_(get_random_key_for_clhash(UINT64_C(0x23a23cf5033c3c81),
                                                     UINT64_C(0xb3816f6a2c68e530)))
    {
    }

    void *get() const {return random_;}

    ~rand_holder()    {free(random_);}
};

const static rand_holder RAND;

namespace std {

  template <typename T>
  struct hash<vector<T>>
  {
    std::uint64_t operator()(const vector<T>& vec) const
    {
        return clhash(RAND.get(), reinterpret_cast<const char *>(vec.data()), vec.size() * sizeof(T));
    }
  };

}


namespace count {

template<typename T>
std::string vec2str(const std::vector<T> &vec) {
    std::string ret;
    for(auto &i: vec) ret += std::to_string(i) + ", ";
    ret.pop_back();
    ret.pop_back();
    return ret;
}

template<typename Test, template<typename...> class Ref>
struct is_specialization : std::false_type {};

template<template<typename...> class Ref, typename... Args>
struct is_specialization<Ref<Args...>, Ref>: std::true_type {};


//template<typename T>

template<typename T, class Hash=std::hash<T>>
class Counter {
    std::size_t                                               n_;
    std::size_t                                           nelem_;
    std::unordered_map<T, std::size_t, Hash>                map_;
    std::unique_ptr<std::unordered_map<unsigned, unsigned>> hist_;

public:
    Counter(): n_(0), nelem_(0), hist_(std::make_unique<std::unordered_map<unsigned, unsigned>>()) {}

    template<typename Q=T>
    void add(const typename std::enable_if<std::is_fundamental<Q>::value, Q>::type elem) {
        auto match(map_.find(elem));
        if(match == map_.end()) map_.emplace(elem, 1);
        else ++match->second;
        ++n_;
    }

    template<typename Q=T>
    void add(const typename std::enable_if<!std::is_fundamental<Q>::value, Q>::type &elem) {
        auto match(map_.find(elem));
        if(match == map_.end()) map_.emplace(elem, 1);
        else ++match->second;
        ++n_;
    }

    void set_nelem(std::size_t nelem) {
        nelem_ = nelem;
    }

    std::size_t get_nelem() const {return nelem_;}

    template<class C>
    void fill(C container) {
        for(auto i: container) add(i);
    }

    std::size_t size()  const {return map_.size();}
    std::size_t total() const {return n_;}
    auto begin()        const {return map_.begin();}
    auto end()          const {return map_.end();}
    auto find(T &elem)  const {return map_.find(elem);}

    const std::unordered_map<T, std::size_t, Hash> &get_map() { return map_;}

    template<typename = std::enable_if<std::is_same<std::vector<std::uint64_t>, T>::value>>
    std::unordered_map<unsigned, unsigned> *make_hist() {
        std::unordered_map<unsigned, unsigned>::iterator m;
        LOG_DEBUG("map size: %zu\n", map_.size());
        for(auto &i: map_) {
            LOG_DEBUG("In map: vec '%s'@%zu\n", vec2str(i.first).data(), i.second);
            if((m = hist_->find(i.second)) == hist_->end()) hist_->emplace(i.second, 1);
            else                                            ++m->second;
        }
        return hist_.get();
    }

    int print_hist(FILE *fp) {
        if(hist_->empty()) make_hist();
        int ret(0);
        std::set<unsigned> countset;
        for(const auto &i: *hist_)
            countset.insert(i.first);
        std::vector<unsigned> counts(countset.begin(), countset.end());
        std::sort(counts.begin(), counts.end());
        fputs("#Count\tNumber of occurrences\n", fp);
        for(auto count: counts) ret += fprintf(fp, "%u\t%u\n", count, hist_->find(count)->second);
        return ret;
    }

    void print_counts(FILE *fp) const {
        struct vecc_t {
            const T *vec_;
            std::size_t count_;
            vecc_t(const T *vec, const std::size_t count): vec_(vec), count_(count) {}
            inline bool operator<(const vecc_t &other) {
                return count_ < other.count_;
            }
        };
        ks::KString ks;
        std::size_t sum(0);
        std::vector<vecc_t> vc;
        vc.reserve(map_.size());
        for(auto &i: map_) vc.emplace_back(&i.first, i.second);
        std::sort(std::begin(vc), std::end(vc));
        ksprintf(ks, "#Number of distinct patterns: %zu\n", vc.size());
        ksprintf(ks, "#Pattern\tCount\tNumber of bits set in pattern\n");
        for(const auto &i: vc) {
            std::size_t n(0);
            const T &vec(*i.vec_);
            std::size_t pc(0);
            for(const auto j: vec) {
                sum += j;
                auto k(j);
                for(std::uint64_t i(0); i < std::min(CHAR_BIT * sizeof(j), nelem_ - n); ++i) {
                    kputc('0' + (k&1), ks);
                    pc += k&1;
                    k >>= 1;
#if 0
                    if(++n >= nelem_) {
                        ksprintf(ks, "\t%zu\n", i.count_);
                        goto end;
                    }
#endif
                }
            }
            ksprintf(ks, "\t%zu\t%zu\n", i.count_, pc);
        }
    
        if(ks.size()) fwrite(ks.data(), 1, ks.size(), fp);
    }
};

} // namespace count

#endif  // COUNTER_H__