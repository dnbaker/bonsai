#ifndef COUNTER_H__
#define COUNTER_H__
#include <algorithm>
#include <typeinfo>
#include <cstddef>
#include <cstdlib>
#include <memory>
#include <set>
#include <unordered_map>
#include <vector>
#include "klib/kstring.h"
#include "lib/logutil.h"
#include "clhash/include/clhash.h"

class rand_holder {
    void *random_;
public:
    rand_holder():
        random_(get_random_key_for_clhash(UINT64_C(0x23a23cf5033c3c81),UINT64_C(0xb3816f6a2c68e530)))
    {
    }
    void *get() const {return random_;}
    ~rand_holder() {free(random_);}
};

const static rand_holder RAND;

namespace std {

  template <typename T>
  struct hash<vector<T>>
  {
    uint64_t operator()(const vector<T>& vec) const
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


template<typename T>
struct vecc_t {
    const T *vec_;
    size_t count_;
    vecc_t(const T *vec, const size_t count): vec_(vec), count_(count) {}
    inline bool operator<(const vecc_t &other) {
        return count_ < other.count_;
    }
};

template<typename T, class Hash=std::hash<T>>
class Counter {
    std::size_t n_;
    std::unordered_map<T, std::size_t, Hash>                map_;
    std::unique_ptr<std::unordered_map<unsigned, unsigned>> hist_;
public:
    Counter(): n_(0), hist_(std::make_unique<std::unordered_map<unsigned, unsigned>>()) {}

    //template<typename = typename std::enable_if<std::is_class<T>::value>::type>
    void add(T &elem) {
        auto match(map_.find(elem));
        if(match == map_.end()) {
            const size_t old(map_.size());
            LOG_DEBUG("New!\n");
            map_.emplace(elem, 1);
            if(map_.size() == old) {
                LOG_DEBUG("Insert failed");
            } else {
                LOG_DEBUG("New size: %zu\n", map_.size());
                for(auto &i: map_) {
                    fprintf(stderr, "first value in key %zu, count: %zu\n", (size_t)i.first[0], i.second);
                }
            }
        } else {
            ++match->second;
            LOG_DEBUG("Found! New count: %zu\n", match->second);
        }
        ++n_;
    }

    template<class C>
    void fill(C container) {
        for(auto i: container) add(i);
    }

    std::size_t size()  const {return map_.size();}
    std::size_t total() const {return n_;}
    auto begin()        const {return map_.begin();}
    auto end()          const {return map_.end();}

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
    //template<typename = std::enable_if<std::is_same<std::vector<std::uint64_t>, T>::value>>
    void print_counts(FILE *fp) {
        fprintf(stderr, "fileno: %i. pointer: %p\n", fileno(fp), (void *)fp);
        kstring_t ks{0, 0, 0};
        size_t sum(0);
        std::vector<vecc_t<T>> vc;
        vc.reserve(map_.size());
        for(auto &i: map_) vc.emplace_back(&i.first, i.second);
        std::sort(std::begin(vc), std::end(vc));
        fprintf(stderr, "Sorted array of vc's of size %zu\n", vc.size());
        for(const auto &i: vc) {
            const std::vector<std::uint64_t> &vec(*i.vec_);
            for(const auto j: vec) {
                sum += j;
                auto k(j);
                for(uint64_t i(0); i < 64; ++i) {
                    kputc('0' + (k&1), &ks);
                    k >>= 1;
                }
            }
            ksprintf(&ks, "\t%zu\n", i.count_);
        }
        if(ks.s) {
            fprintf(stderr, "Made buf: %s\n", ks.s);
            fwrite(ks.s, 1, ks.l, fp);
            free(ks.s);
        }
        fprintf(stderr, "Max count: %zu. Min: %zu\n", vc.size() ? vc[vc.size() - 1].count_: 0ul, vc.size() ? vc[0].count_: 0ul);
    }
};


template<typename T>
size_t unique(T &vec, Counter<uint64_t> &c);

} // namespace count

#endif  // COUNTER_H__
