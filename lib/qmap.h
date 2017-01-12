#ifndef QMAP_H_
#define QMAP_H_

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
    T el_;
    ScoreType score_;
    INLINE ElScore(T el, ScoreType score): el_(el), score_(score) {
    }
    INLINE ElScore(): el_(0), score_(0) {}
    INLINE bool operator <(const ElScore &other) const {
        return std::tie(score_, el_) < std::tie(other.score_, other.el_);
    }
    INLINE bool operator ==(const ElScore &other) const {
        return el_ == other.el_;
        //return std::tie(score_, el_) == std::tie(other.el_, other.score_); // Lexicographic is tie-breaker.
    }
};


template<typename T, typename ScoreType>
class QueueMap {
    // I could make this more efficient by using pointers instead of
    // ElScore structs.
    // *maybe* TODO
    // Could also easily templatify this module for other windowing tasks.
    typedef typename std::map<ElScore<T, ScoreType>, unsigned>::iterator map_iterator;
    std::list<ElScore<T, ScoreType>> list_;
#if !NDEBUG
public:
    std::map<ElScore<T, ScoreType>, unsigned> map_;
private:
#else
    std::map<ElScore<T, ScoreType>, unsigned> map_;
#endif
    const size_t wsz_;  // window size to keep
    public:
    INLINE void add(const ElScore<T, ScoreType> &el) {
        auto it(map_.lower_bound(el));
        if(it != map_.end()) {
            if(it->first == el) ++it->second;
            else map_.emplace_hint(it, el, 1);
        }
        else map_.emplace(el, 1);
        //else map_.emplace_hint(it, el, 1);
    }
    INLINE void del(const ElScore<T, ScoreType> &el) {
        auto f(map_.find(el));
        //LOG_DEBUG("Removing %s\n", f->first.to_string().data());
        if(--f->second <= 0)
            map_.erase(f);
    }
    map_iterator begin() {
        return map_.begin();
    }
    map_iterator end() {
        return map_.end();
    }
    uint64_t next_value(const T el, const uint64_t score) {
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
    QueueMap(size_t wsz): wsz_(wsz) {}
    void reset() {
        list_.clear();
        map_.clear();
    }
};

using qmap_t = QueueMap<uint64_t, uint64_t>;
using elscore_t = ElScore<uint64_t, uint64_t>;

} // namespace emp

#endif //ifndef QMAP_H_
