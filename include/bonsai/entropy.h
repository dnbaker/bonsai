#pragma once
#include "circularqueue/cq.h"
#include "bonsai/kmerutil.h"
#include "flat_hash_map/flat_hash_map.hpp"

namespace bns {

class CircusEnt {
    using CountT = uint64_t;
    circ::deque<char> q_;
    size_t        cqsz_;
    const size_t   qsz_;
    const double qszinv_;
    ska::flat_hash_map<char, size_t> counts_;
public:
    static constexpr double NOT_FULL = -1.;
    CircusEnt(size_t qsz): q_(qsz), cqsz_(0), qsz_(qsz), qszinv_(1./qsz) {
    }
    CircusEnt(const CircusEnt &other): q_(other.q_), cqsz_(other.cqsz_), qsz_(other.qsz_), qszinv_(other.qszinv_), counts_(other.counts_) {
    }
    CircusEnt(CircusEnt &&other) = default;
    void clear() {
        counts_.clear();
        cqsz_ = 0;
        q_.clear();
    }
    void push(char c) {
        // Increment count of added nucleotide
        // If not full
        auto it = counts_.find(c);
        if(it == counts_.end()) it = counts_.emplace(c, 1).first;
        else ++it->second;
        if(cqsz_ == qsz_) {
            it = counts_.find(q_.pop());
            assert(it != counts_.end());
            if(--it->second == 0)
                counts_.erase(it);
            // Pop and decrement
        } else
            ++cqsz_; // Or just keep filling the window
        q_.push(c);
    }
    double value() const {
        if(unlikely(cqsz_ < qsz_)) return NOT_FULL;
        return std::accumulate(counts_.begin(), counts_.end(), 0.,
               [qi=qszinv_](double s, auto v) {return s + v.second * qi * std::log(v.second * qi);});
    }
    double next_ent(char c) {
        push(c);
        return value();
    }
};

} // namespace bns
