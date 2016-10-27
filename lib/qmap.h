#include <cstdint>
#include <cstdio>
#include <cinttypes>
#include <map>
#include <list>
#include <vector>

#include "util.h"

namespace kpg {

struct elscore_t {
    uint64_t el_, score_;
    elscore_t(uint64_t el, uint64_t score): el_(el), score_(score) {}
    elscore_t(): el_(0), score_(0) {}
    inline bool operator <(const elscore_t &other) const {
        return score_ < other.score_ || el_ < other.el_; // Lexicographic is tie-breaker.
    }
};

struct esq_t: public std::list<elscore_t> {
};

struct esmap_t: public std::map<elscore_t, unsigned> {
    void add(const elscore_t &el) {
        auto f(find(el));
        if(f == end())
            emplace(el, 1);
        else
            ++f->second;
    }
    void del(const elscore_t &el) {
        const auto f(find(el));
        if(f != end())
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
    esmap_t map_;
    const size_t wsz_;  // window size to keep
    public:
    uint64_t next_value(const uint64_t el, const uint64_t score) {
        if(list_.size() == wsz_) {
            map_.del(*list_.begin());
            list_.pop_front();
        }
        list_.emplace_back(el, score);
        map_.add(*list_.end());
        return list_.size() == wsz_ ? map_.begin()->first.el_: (uint64_t)~0;
        // Signal a window that is not filled by 0xFFFFFFFFFFFFFFFF
    }
    qmap_t(size_t wsz): wsz_(wsz) {}
    void reset() {
        list_.clear();
        map_.clear();
    }
};

}
