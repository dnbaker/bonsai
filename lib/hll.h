#ifndef _HLL_H_
#define _HLL_H_
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <cinttypes>
#include <cassert>
#include <algorithm>
#include <vector>

#include "hash.h"

namespace kpg {

constexpr double make_alpha(size_t m) {
    switch(m) {
        case 16: return .673;
        case 32: return .697;
        case 64: return .709;
        default: return 0.7213 / (1 + 1.079/m);
    }
}

template<size_t np>
class hll_t {
// This is a HyperLogLog implementation.
// To make it general, the actual point of entry is a 64-bit integer hash function.
// Therefore, you have to perform a hash function to convert various types into a suitable query.
// I've added an overload for strings using others' hash functions for convenience.

// Attributes
public:
    static const size_t m_ = 1 << np;
    static constexpr double alpha_ = make_alpha(m_);
    static constexpr double relative_error_ = 1.03896 / std::sqrt(m_);
    static const uint64_t bitmask_ = (~((uint64_t)0)) >> np;
// Members
    std::vector<uint8_t> core_;
private:
    double sum_;
    int is_calculated_;
// Functions
public:
    // Call sum to recalculate if you have changed contents.
    void sum() {
        sum_ = 0;
        for(unsigned i(0); i < m_; ++i) sum_ += 1. / (1 << core_[i]);
        is_calculated_ = 1;
    }
    double report() {
        if(!is_calculated_) sum();
        const double ret = alpha_ * m_ * m_ / sum_;
        // Correct for small values
        if(ret < m_ * 2.5) {
            int t(0);
            for(int i(0); i < m_; ++i) t += (core_[i] == 0);
            if(t) return m_ * std::log((double)(m_) / t);
        }
        return ret;
        // We don't correct for too large just yet, but we should soon.
    }
    double est_err() {
        return relative_error_ * report();
    }
    void add(uint64_t hashval) {
        const uint32_t index = hashval >> (64 - np);
        const int lzt(__builtin_clzll(hashval << np) + 1);
        assert(index < m_);
        assert(!hashval || lzt == __builtin_clzll(hashval << np) + 1);
        if(core_[index] < lzt) core_[index] = lzt;
    }
    void add(const char *str) {
        add(((uint64_t)X31_hash_string(str) << 32) | dbm_hash(str));
    }
    void add(const char *str, size_t l) {
        add(((uint64_t)X31_hash_string(str) << 32) | dbm_hash(str, l));
    }
    hll_t(): core_(m_, 0), sum_(0.), is_calculated_(0) {}
    ~hll_t() {}
    hll_t<np> const &operator+=(const hll_t<np> &other) {
         // If we ever find this to be expensive, this could be trivially implemented with SIMD.
        for(unsigned i(0); i < m_; ++i) core_[i] |= other.core_[i];
        return *this;
    }
    std::string to_string() {
        return std::to_string(report()) + ", +- " + std::to_string(est_err());
    }
    // Note: We store values as (64 - value) and take the maximum when updating,
    // which allows us to minimize our trailing zeroes without sacrificing calloc.
};

typedef hll_t<22> HyperLogLog;

} // namespace kpg

#endif // #ifndef _HLL_H_
