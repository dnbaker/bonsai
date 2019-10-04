#pragma once
#include "cq.h"
#include "kmerutil.h"

namespace bns {

class CircusEnt {
    circ::deque<char> q_;
    uint8_t   counts_[4];
    uint8_t        cqsz_;
    const uint8_t   qsz_;
    const double qszinv_;
public:
    static constexpr double NOT_FULL = -1.;
    CircusEnt(unsigned qsz): q_(qsz), counts_{0}, cqsz_(0), qsz_(qsz), qszinv_(1./qsz) {
        if(qsz == 0 || qsz > 255) RUNTIME_ERROR(std::string("Illegal queue size with qsz = " + std::to_string(qsz) + " and maximum value " + std::to_string(std::numeric_limits<uint8_t>::max())));
    }
    CircusEnt(const CircusEnt &other): q_(other.q_), qsz_(other.qsz_), qszinv_(other.qszinv_) {
        std::memcpy(counts_, other.counts_, sizeof(counts_));
    }
    CircusEnt(CircusEnt &&other) = default;
    void clear() {
        *(reinterpret_cast<uint32_t *>(counts_)) = cqsz_ = 0;
        q_.clear();
    }
    void push(char c) {
#if !NDEBUG
        if(cqsz_ == qsz_) {
            auto to_pop = cstr_lut[q_.pop()];
            assert(to_pop != static_cast<int8_t>(-1));
            assert(counts_[to_pop]);
            --counts_[to_pop];
        }
        q_.push(c);
        assert(cstr_lut[c] != static_cast<int8_t>(-1));
        ++counts_[cstr_lut[c]];
        static_assert(std::is_unsigned<typename std::decay<decltype(q_.size())>::type>::value, "q's size should be unsigned.");
        assert(unsigned(counts_[0] + counts_[1] + counts_[2] + counts_[3]) == q_.size());
#else
        // Increment count of added nucleotide
        // If not fu
        ++counts_[cstr_lut[c]];
        if(cqsz_ == qsz_) {
            --counts_[cstr_lut[q_.pop()]]; // Pop and decrement
        } else {
            ++cqsz_;                       // Or just keep filling the window
        }
        q_.push(c);
#endif
    }
    double value() const {
        if(unlikely(cqsz_ < qsz_)) return NOT_FULL;
        assert(cqsz_ == q_.size());
        double tmp(qszinv_ * (counts_[0])), sum(tmp * std::log2(tmp));
        tmp = qszinv_ * counts_[1], sum += tmp * std::log2(tmp);
        tmp = qszinv_ * counts_[2], sum += tmp * std::log2(tmp);
        tmp = qszinv_ * counts_[3], sum += tmp * std::log2(tmp);
        return sum;
    }
    u64 score() const {
        return UINT64_MAX - static_cast<u64>(UINT64_C(7958933093282078720) * value());
    }
    double next_ent(char c) {
        push(c);
        return value();
    }
};

} // namespace bns
