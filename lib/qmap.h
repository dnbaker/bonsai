#ifndef _QMAP_H__
#define _QMAP_H__

#include <cstdint>
#include <cstdio>
#include <cinttypes>
#include <map>
#include <list>
#include <vector>

#include "util.h"
#include "kmerutil.h"

namespace emp {

template<typename T, typename ScoreType>
struct ElScore {
    T            el_;
    ScoreType score_;
    INLINE ElScore(T el, ScoreType score): el_(el), score_(score) {
    }
    INLINE ElScore(): el_(0), score_(0) {}
    INLINE bool operator<(const ElScore &other) const {
        return std::tie(score_, el_) < std::tie(other.score_, other.el_);
    }
    INLINE bool operator==(const ElScore &other) const {
        return el_ == other.el_;
        //return std::tie(score_, el_) == std::tie(other.el_, other.score_); // Lexicographic is tie-breaker.
    }
};


template<typename T, typename ScoreType>
class QueueMap {
    // I could make this more efficient by using pointers instead of
    // ElScore structs.
    // *maybe* TODO
    //
    using PairType = ElScore<T, ScoreType>;
    using map_iterator = typename std::map<ElScore<T, ScoreType>, unsigned>::iterator;
    std::list<ElScore<T, ScoreType>>         list_;
    std::map<ElScore<T, ScoreType>, unsigned> map_;
    const std::size_t                         wsz_;  // window size to keep
    public:
    QueueMap(std::size_t wsz): wsz_(wsz) {}
    INLINE void add(const PairType &el) {
        auto it(map_.lower_bound(el));
        if(it != map_.end()) {
            if(it->first == el) ++it->second;
            else map_.emplace_hint(it, el, 1);
        }
        else map_.emplace(el, 1);
        //else map_.emplace_hint(it, el, 1);
    }
    INLINE void del(const PairType &el) {
        auto f(map_.find(el));
        //LOG_DEBUG("Removing %s\n", f->first.to_string().data());
        if(--f->second <= 0)
            map_.erase(f);
    }
    map_iterator begin() {
        return map_.begin();
    }
    const map_iterator begin() const {
        return map_.cbegin();
    }
    map_iterator end() {
        return map_.end();
    }
    const map_iterator end() const {
        return map_.cend();
    }
    // Do a std::enable_if that involves moving the element if it's by reference?
    std::uint64_t next_value(const T el, const std::uint64_t score) {
        list_.emplace_back(el, score);
        add(list_.back());
        if(list_.size() > wsz_) {
            //fprintf(stderr, "list size: %zu. wsz: %zu\n", list_.size(), wsz_);
            //map_.del(list_.front());
            del(list_.front());
            list_.pop_front();
        }
        return list_.size() == wsz_ ? map_.begin()->first.el_: BF;
        // Signal a window that is not filled by 0xFFFFFFFFFFFFFFFF
    }
    void reset() {
        list_.clear();
        map_.clear();
    }
};

using qmap_t = QueueMap<std::uint64_t, std::uint64_t>;
using elscore_t = ElScore<std::uint64_t, std::uint64_t>;

} // namespace emp

#endif //ifndef _QMAP_H__
