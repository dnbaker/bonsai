#ifndef _QMAP_H__
#define _QMAP_H__

#include <cstdint>
#include <cstdio>
#include <cinttypes>
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
        return std::tie(score_, el_) < std::tie(other.score_, el_); // Lexicographic is tie-breaker.
    }
    INLINE bool operator==(const ElScore &other) const {
        return el_ == other.el_;
    }
};


/*
Sliding-window minimum over (score, element) pairs with exact first-in
tie-breaking. The queue is kept monotonically increasing from front to
back: a new pair pops every back pair that it does not strictly exceed,
so the front always holds the lexicographic minimum of the window —
identical selection to the previous btree-multimap implementation
(minimum score, smallest element on ties) at amortized O(1) per slide
instead of O(log w) tree operations per insert and erase. Window fullness
is tracked by position: dominance pops shrink the deque below the window
population, so expiry and fullness use sequence numbers, not deque size.
*/
template<typename T, typename ScoreType, typename Compare=std::less<void>>
class QueueMap {
    using PairType           = ElScore<T, ScoreType>;
    using list_type = circ::deque<ElScore<T, ScoreType>, u32>;
    using pos_type  = circ::deque<u64, u32>;

    list_type list_;  // monotone: front = min (score, el) of live window entries
    pos_type  poss_;  // sequence number of each list_ entry, aligned 1:1
    u64 seq_;         // number of pushes since construction/reset
    size_t wsz_;      // window size to keep
    public:
    QueueMap(size_t wsz=1): list_(wsz), poss_(wsz), seq_(0), wsz_(wsz) {}
    QueueMap(QueueMap &&o) = default;
    QueueMap &operator=(QueueMap &&o) = default;
    QueueMap &operator=(const QueueMap &o) = default;
    void resize(size_t newsz) {
        wsz_ = newsz;
        list_.resize(newsz);
        poss_.resize(newsz);
    }
    // Monotone push: back pairs not strictly smaller than el are dead —
    // el outscores them and expires later.
    INLINE void add(const PairType &el) {
        while(list_.size() && !(list_.back() < el)) {
            list_.pop_back();
            poss_.pop_back();
        }
        list_.emplace_back(el);
        poss_.emplace_back(seq_++);
    }
    // The only expiry pattern is the window sliding past the front pair.
    INLINE void del(const PairType &) {
        list_.pop_front();
        poss_.pop_front();
    }
    INLINE T next_value(const T el, const T score) {
        add(PairType(el, score));
        // Expire the element wsz positions back; if dominance already removed
        // it, the front is younger and stays.
        if(seq_ > wsz_ && poss_.front() == seq_ - wsz_ - 1) {
            list_.pop_front();
            poss_.pop_front();
        }
        if(seq_ >= wsz_) return list_.front().el_;
        if(std::is_same<T, u128>::value)
            return u128(-1);
        return std::numeric_limits<T>::max();
        // Signal a window that is not filled by 0xFFFFFFFFFFFFFFFF
    }
    void reset() {
        list_.clear();
        poss_.clear();
        seq_ = 0;
    }
    size_t size() const {return wsz_;}
    size_t n_in_queue() const {return std::min(seq_, u64(wsz_));}
    bool partially_full() const {
        const auto n = n_in_queue();
        return n > 0u && n < wsz_;
    }
    // The minimum-(score, el) pair currently in the window.
    const auto &max_in_queue() const {return list_.front();}
};

using qmap_t = QueueMap<u64, u64>;
using qmapf_t = QueueMap<u64, double>;
using qmap128_t = QueueMap<u128, u64>;
using qmap128f_t = QueueMap<u128, double>;
using elscore_t = ElScore<u64, u64>;

} // namespace bns

#endif //ifndef _QMAP_H__
