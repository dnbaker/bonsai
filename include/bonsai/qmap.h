#ifndef _QMAP_H__
#define _QMAP_H__

#include <cstdint>
#include <cstdio>
#include <cinttypes>
#include "btree/map.h"
#include "hll/include/circularqueue/cq.h"
#include <vector>

#include "bonsai/util.h"
#include "bonsai/kmerutil.h"

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


template<typename T, typename ScoreType, typename Compare=std::less<void>>
class QueueMap {
    using PairType           = ElScore<T, ScoreType>;
    using MapType = btree::map<PairType, unsigned>;
    using map_iterator = typename MapType::iterator;
    using const_map_iterator = typename MapType::const_iterator;
    using list_type = circ::deque<ElScore<T, ScoreType>, u32>;

    list_type list_;
    MapType map_;
    size_t wsz_;  // window size to keep
    public:
    QueueMap(size_t wsz=1): list_(wsz), wsz_(wsz) {}
    QueueMap(QueueMap &&o): list_(std::move(o.list_)), map_(std::move(o.map_)), wsz_(o.wsz_) {}
    QueueMap(const QueueMap &&o): list_(o.list_), map_(o.map_), wsz_(o.wsz_) {}
    QueueMap &operator=(QueueMap &&o) {
        list_ = std::move(o.list_);
        map_ = std::move(o.map_);
        wsz_ = std::move(o.wsz_);
        return *this;
    }
    QueueMap &operator=(const QueueMap &o) {
        list_ = o.list_;
        map_ = o.map_;
        wsz_ = o.wsz_;
        return *this;
    }
    void resize(size_t newsz) {
        wsz_ = newsz;
        list_.resize(newsz);
    }
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
    T next_value(const T el, const T score) {
        add(list_.emplace_back(el, score));
        if(list_.size() > wsz_) del(list_.pop_front());
        if(list_.size() == wsz_) return map_.begin()->first.el_;
        if(std::is_same<T, u128>::value)
            return u128(-1);
        return std::numeric_limits<T>::max();
        // Signal a window that is not filled by 0xFFFFFFFFFFFFFFFF
    }
    void reset() {
        list_.clear(); map_.clear();
    }
    size_t size() const {return wsz_;}
    size_t n_in_queue() const {return list_.size();}
    bool partially_full() const {
        const auto n = n_in_queue();
        return n > 0u && n < wsz_;
    }
    const auto &max_in_queue() const {return begin()->first;}
};

using qmap_t = QueueMap<u64, u64>;
using qmapf_t = QueueMap<u64, double>;
using qmap128_t = QueueMap<u128, u64>;
using qmap128f_t = QueueMap<u128, double>;
using elscore_t = ElScore<u64, u64>;

} // namespace bns

#endif //ifndef _QMAP_H__
