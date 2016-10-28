#include <cstdint>
#include <cstdio>
#include <cinttypes>
#include <map>
#include <list>
#include <vector>

#include "util.h"
#include "kmerutil.h"

namespace kpg {

struct elscore_t {
    uint64_t el_, score_;
    INLINE elscore_t(uint64_t el, uint64_t score): el_(el), score_(score) {
        assert(el == el_);
        assert(score == score_);
    }
    INLINE elscore_t(): el_(0), score_(0) {}
    inline bool operator <(const elscore_t &other) const {
        return score_ < other.score_ || el_ < other.el_; // Lexicographic is tie-breaker.
    }
    inline bool operator ==(const elscore_t &other) const {
        return score_ == other.score_ && el_ == other.el_; // Lexicographic is tie-breaker.
    }
};

struct esq_t: public std::list<elscore_t> {
};

struct esmap_t: public std::map<elscore_t, unsigned> {
    std::pair<std::map<elscore_t, unsigned>::iterator, bool> finder_;
    void add(const elscore_t &el) {
        std::map<elscore_t, unsigned>::iterator it(lower_bound(el));
        if(it->first == el) ++it->second;
        else emplace_hint(it, el, 1);
    }
    void del(const elscore_t &el) {
        const auto f(find(el));
        assert(f != end());
        if(--f->second <= 0)
            erase(f);
    }
};

class qmap_t {
    // I could make this more efficient by using pointers instead of
    // elscore_t structs.
    // *maybe* TODO
    // Could also easily templatify this module for other windowing tasks.
    esq_t list_;
#if !NDEBUG
public:
    esmap_t map_;
private:
#endif
    const size_t wsz_;  // window size to keep
    public:
    uint64_t next_value(const uint64_t el, const uint64_t score) {
        list_.emplace_back(el, score);
        assert(list_.back().el_ == el);
        assert(list_.back().score_ == score);
        map_.add(list_.back());
        if(list_.size() > wsz_) {
            map_.del(list_.front());
            list_.pop_front();
        }
        assert(list_.size() <= wsz_);
        return list_.size() == wsz_ ? map_.begin()->first.el_: BF;
        // Signal a window that is not filled by 0xFFFFFFFFFFFFFFFF
    }
    qmap_t(size_t wsz): wsz_(wsz) {}
    void reset() {
        list_.clear();
        map_.clear();
    }
};

}
