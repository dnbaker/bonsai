#ifndef _COUNTER_H__
#define _COUNTER_H__
#include <algorithm>
#include <typeinfo>
#include <cstddef>
#include <cstdlib>
#include <climits>
#include <memory>
#include <set>
#include <unordered_map>
#include <vector>
#include "util.h"
#include "kspp/ks.h"
#include "hash.h"

namespace emp {

class rand_holder {
    const void *random_;
public:
    rand_holder(): random_(get_random_key_for_clhash(UINT64_C(0x23a23cf5033c3c81),
                                                     UINT64_C(0xb3816f6a2c68e530)))
    {
    }

    const void *get()        const {return random_;}
    const void *operator->() const {return random_;}
    operator const void *()  const {return random_;}
    template<typename T>
    u64 hash(const T *p, size_t nbytes) {
        return clhash(random_, p, nbytes);
    }
    template<typename T>
    u64 hash(const std::vector<T>& vec) {
        return hash(vec.data(), vec.size() * sizeof(T));
    }

    ~rand_holder()    {std::free(const_cast<void *>(random_));}
};

const static rand_holder RAND;

} // namespace emp

namespace std {

  template <typename T, typename size_type>
  struct hash<lazy::vector<T, size_type>>
  {
    uint64_t operator()(const lazy::vector<T, size_type>& vec) const
    {
        return clhash(emp::RAND, reinterpret_cast<const char *>(vec.data()), vec.size() * sizeof(T));
    }
  };
  template <typename T>
  struct hash<vector<T>>
  {
    uint64_t operator()(const vector<T>& vec) const
    {
        return clhash(emp::RAND, reinterpret_cast<const char *>(vec.data()), vec.size() * sizeof(T));
    }
  };


  template <typename T, typename size_type>
  struct hash<pair<int, lazy::vector<T, size_type>>>
  {
    uint64_t operator()(const pair<int, lazy::vector<T, size_type>>& p) const
    {
        return clhash(emp::RAND, reinterpret_cast<const char *>(p.second.data()), p.second.size() * sizeof(T)) ^ ((u64(p.first) << 32) | p.first);
    }
  };

  template <typename T>
  struct hash<pair<int, vector<T>>>
  {
    uint64_t operator()(const pair<int, vector<T>>& p) const
    {
        return clhash(emp::RAND, reinterpret_cast<const char *>(p.second.data()), p.second.size() * sizeof(T)) ^ ((u64(p.first) << 32) | p.first);
    }
  };

}

namespace emp {

namespace count {

template<typename Container, typename=std::enable_if_t<std::is_arithmetic_v<typename Container::value_type>>>
std::string vec2str(const Container &vec)  {
    std::string ret;
    for(const auto i: vec) ret += std::to_string(i) + ", ";
    ret.pop_back(), ret.pop_back();
    return ret;
}

template<typename Test, template<typename...> class Ref>
struct is_specialization : std::false_type {};

template<template<typename...> class Ref, typename... Args>
struct is_specialization<Ref<Args...>, Ref>: std::true_type {};

template<typename T, class Hash=std::hash<T>, typename SizeType=size_t, typename=std::enable_if_t<std::is_arithmetic_v<SizeType>>>
class Counter {
    SizeType                                     n_;
    SizeType                                 nelem_;
    std::unordered_map<T, SizeType, Hash>      map_;
    std::unordered_map<unsigned, unsigned>   *hist_;

public:
    Counter(): n_(0), nelem_(0), hist_(nullptr) {}
    ~Counter() {if(hist_) delete hist_;}

    template<typename Q=T>
    void add(const typename std::enable_if_t<std::is_fundamental_v<Q>, Q> elem) {
#if __GNUC__ >= 7
        if(auto match = map_.find(elem); match == map_.end()) map_.emplace(elem, 1);
        else                                          ++match->second;
#else
        auto match(map_.find(elem));
        if(match == map_.end()) map_.emplace(elem, 1);
        else                    ++match->second;
#endif
        ++n_;
    }

    template<typename Q=T>
    void add(const typename std::enable_if_t<!std::is_fundamental_v<Q>, Q> &elem) {
        auto match(map_.find(elem));
        if(match == map_.end()) map_.emplace(elem, 1);
        else                    ++match->second;
        ++n_;
    }
    template<typename Q=T>
    void add(typename std::enable_if_t<!std::is_fundamental_v<Q>, Q> &&elem) {
        auto match(map_.find(elem));
        if(match == map_.end()) map_.emplace(std::move(elem), 1);
        else                    ++match->second;
        ++n_;
    }

    void set_nelem(SizeType nelem) {
        nelem_ = nelem;
    }

    SizeType get_nelem() const {return nelem_;}

    template<class C>
    void fill(C container) {
        for(const auto i: container) add(i);
    }
    template<class C, class Fn>
    void fill(C container, const Fn &function) {
        for(const auto i: container) add(function(i));
    }

    SizeType size()        const {return map_.size();}
    SizeType total()       const {return n_;}
    auto begin()              const {return map_.begin();}
    auto begin()                    {return map_.begin();}
    auto end()                const {return map_.end();}
    auto end()                      {return map_.end();}
    auto cbegin()             const {return map_.cbegin();}
    auto cend()               const {return map_.cend();}
    auto find(const T &elem)  const {return map_.find(elem);}

    const std::unordered_map<T, SizeType, Hash> &get_map() const { return map_;}

    template<typename = std::enable_if<std::is_same_v<bitvec_t, T>>>
    std::unordered_map<unsigned, unsigned> *make_hist() {
        if(hist_) delete hist_;
        hist_ = new std::unordered_map<unsigned, unsigned>;
        std::unordered_map<unsigned, unsigned>::iterator m;
        LOG_DEBUG("map size: %zu\n", map_.size());
        for(auto &i: map_) {
            LOG_DEBUG("In map: vec '%s'@%zu\n", vec2str(i.first).data(), i.second);
            if((m = hist_->find(i.second)) == hist_->end()) hist_->emplace(i.second, 1);
            else                                            ++m->second;
        }
        return hist_;
    }

    int print_hist(std::FILE *fp) {
        if(hist_ == nullptr || hist_->empty()) make_hist();
        int ret(0);
        std::set<unsigned> countset;
        for(const auto &i: *hist_)
            countset.insert(i.first);
        std::vector<unsigned> counts(countset.begin(), countset.end());
        SORT(counts.begin(), counts.end(), [](auto a, auto b) {return a < b;});
        std::fputs("#Count\tNumber of occurrences\n", fp);
        for(auto count: counts) ret += std::fprintf(fp, "%u\t%u\n", count, hist_->find(count)->second);
        return ret;
    }

    void print_vec() const {
        for(const auto &i: map_) {
            std::fprintf(stderr, "%s@%zu\n", vec2str(i.first).data(), i.second);
        }
    }

    void print_counts(std::FILE *fp) const {
        if constexpr(!std::is_scalar_v<T>) {
            struct vecc_t {
                const T *vec_;
                SizeType count_;
                vecc_t(const T *vec, const SizeType count): vec_(vec), count_(count) {}
#if 0
                inline bool operator<(const vecc_t &other) {
                    return count_ < other.count_;
                }
#endif
            };
            ks::string ks;
            ks.resize(256);
            SizeType sum(0);
            std::vector<vecc_t> vc;
            vc.reserve(map_.size());
            LOG_DEBUG("Size of map: %zu\n", map_.size());
            for(const auto &i: map_) {
                vc.emplace_back(&i.first, i.second);
            }
            pdqsort(std::begin(vc), std::end(vc), [](const auto &x, const auto &y) {return x.count_ < y.count_;});
            ks.sprintf("#Number of distinct patterns: %zu\n", vc.size());
            ks.sprintf("#Pattern\tCount\tNumber of bits set in pattern\n");
            if(map_.empty()) return;
            for(const auto &i: vc) {
                SizeType n(0), pc(0);
                std::fprintf(stderr, "Pointer: %p\n", (void *)i.vec_);
                std::fprintf(stderr, "Size: %zu\n", i.vec_->size());
                for(auto j: *i.vec_) {
                    std::fprintf(stderr, "Element %" PRIu64 " in vector.\n", j);
                    sum += j;
                    for(u64 i(0); i < std::min(CHAR_BIT * sizeof(j), nelem_ - n); ++i)
                        ks.putc_('0' + (j&1)), pc += j&1, j >>= 1;
                }
                ks.sprintf("\t%zu\t%zu\n", i.count_, pc);
            }
            if(ks.size()) std::fwrite(ks.data(), 1, ks.size(), fp);
        } else {
            for(const auto &el: map_) {
                std::fprintf(fp, "%s\t%s\n", std::to_string(el.first).data(), std::to_string(el.second).data());
            }
        }
    }
};

} // namespace count

} // namespace emp

#endif  // _COUNTER_H__
