#ifndef COUNTER_H__
#define COUNTER_H__
#include <unordered_map>
#include <cstddef>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include "klib/kstring.h"

namespace count {

template<typename T, class Hash=std::hash<T>>
class Counter {
    std::unordered_map<T, std::size_t, Hash> map_;
    std::size_t n_;
public:
    Counter(): n_(0) {}

#if 0
    void add(T elem) {
        auto match(map_.find(elem));
        if(match == map_.end()) map_.emplace(elem, 1);
        else ++match->second;
        ++n_;
    }
#endif

    //template<typename = typename std::enable_if<std::is_class<T>::value>::type>
    void add(T &elem) {
        auto match(map_.find(elem));
        if(match == map_.end()) map_.emplace(elem, 1);
        else ++match->second;
        ++n_;
    }
#if 0
    template<typename = typename std::enable_if<!std::is_class<T>::value>::type>
    void add(T elem) {
        auto match(map_.find(elem));
        if(match == map_.end()) map_.emplace(elem, 1);
        else ++match->second;
        ++n_;
    }
#endif

    template<class C>
    void fill(C container) {
        for(auto i: container) add(i);
    }
    std::size_t size()  const {return map_.size();}
    std::size_t total() const {return n_;}
    auto begin()        const {return map_.begin();}
    auto end()          const {return map_.end();}
    template<typename = std::enable_if<std::is_same<std::vector<std::uint64_t>, T>::value>>
    void print_counts(FILE *fp) {
        kstring_t ks{0, 0, 0};
        size_t sum(0);
        struct vecc_t {
            const std::vector<std::uint64_t> *vec_;
            size_t count_;
            vecc_t(const std::vector<std::uint64_t> &vec, size_t count): vec_(&vec), count_(count) {}
        };
        std::vector<vecc_t> vc;
        vc.reserve(map_.size());
        for(auto &i: map_) vc.emplace_back(i.first, i.second);
        std::sort(std::begin(vc), std::end(vc), [](vecc_t &a, vecc_t &b) {
            return a.count_ < b.count_;
        });
        fprintf(stderr, "Max count: %zu. Min: %zu\n", vc[vc.size() - 1].count_, vc[0].count_);
        for(auto &i: vc) {
            for(auto j: *i.vec_) {
                sum += j;
                auto k(j);
                for(uint64_t i(0); i < 64; ++i) {
                    kputc('0' + (k&1), &ks);
                    k >>= 1;
                }
            }
            ksprintf(&ks, "\t%zu\n", i.count_);
        }
        fwrite(ks.s, 1, ks.l, fp);
        free(ks.s);
        //fprintf(stderr, "Sum: %zu\n", sum);
    }
};

template<typename T>
size_t unique(T &vec, Counter<uint64_t> &c);

} // namespace count

#endif  // COUNTER_H__
