#ifndef _QMAP_H__
#define _QMAP_H__

#include <cstdint>
#include <cstdio>
#include <cinttypes>
#include <map>
#include "circular_buffer.h"
#include <vector>

#include "util.h"
#include "kmerutil.h"

namespace bns {

template<typename T, typename ScoreType>
struct ElScore {
    T            el_;
    ScoreType score_;
    INLINE ElScore(T el, ScoreType score): el_(el), score_(score) {
    }
    INLINE ElScore(): el_(0), score_(0) {}
    INLINE bool operator<(const ElScore &other) const {
        return std::tie(score_, el_) < std::tie(other.score_, other.el_); // Lexicographic is tie-breaker.
    }
    INLINE bool operator==(const ElScore &other) const {
        return el_ == other.el_;
    }
};


template<typename T, typename ScoreType>
class QueueMap {
    // I could make this more efficient by using pointers instead of
    // ElScore structs.
    //
    using PairType           = ElScore<T, ScoreType>;
    using map_iterator       = typename std::map<ElScore<T, ScoreType>, unsigned>::iterator;
    using const_map_iterator = typename std::map<ElScore<T, ScoreType>, unsigned>::const_iterator;
    using list_type = circ::deque<ElScore<T, ScoreType>, u32>;

    list_type list_;
    std::map<ElScore<T, ScoreType>, u32> map_;
    const size_t                         wsz_;  // window size to keep
    public:
    QueueMap(size_t wsz): list_(wsz), wsz_(wsz) {}
    INLINE void add(const PairType &el) {
        auto it(map_.lower_bound(el));
        if(it != map_.end()) {
            if(it->first == el) ++it->second;
            else map_.emplace_hint(it, el, 1);
        } else map_.emplace(el, 1);
    }
    INLINE void del(const PairType &el) {
        auto f(map_.find(el));
        if(--f->second <= 0) map_.erase(f);
    }
    map_iterator       begin()       {return map_.begin();}
    const_map_iterator begin() const {return map_.cbegin();}
    map_iterator       end()       {return map_.end();}
    const_map_iterator end() const {return map_.cend();}
    // Do a std::enable_if that involves moving the element if it's by reference?
    u64 next_value(const T el, const u64 score) {
        add(list_.emplace_back(el, score));
        if(list_.size() > wsz_) del(list_.pop_front());
        return list_.size() == wsz_ ? map_.begin()->first.el_: BF;
        // Signal a window that is not filled by 0xFFFFFFFFFFFFFFFF
    }
    void reset() {
        list_.clear();
        map_.clear();
    }
};

using qmap_t = QueueMap<u64, u64>;
using elscore_t = ElScore<u64, u64>;

} // namespace bns

#endif //ifndef _QMAP_H__
