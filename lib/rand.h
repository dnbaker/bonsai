#ifndef _EMP_RAND_H__
#define _EMP_RAND_H__
#include <random>
#include "fast_mutex.h"

namespace rng {

struct RandTwister {
    std::mt19937_64 twister_;

    using ResultType     = std::mt19937_64::result_type;

    static const ResultType MAX     = std::mt19937_64::max();
    static const ResultType MIN     = std::mt19937_64::min();
    static constexpr double MAX_INV = 1. / MAX;

    RandTwister(ResultType seed=std::rand()): twister_(seed) {}
    auto operator()()                              {return twister_();}
    auto operator()(std::mt19937_64 &engine) const {return engine();}
    // Generate a large number of random integers.
    void operator()(size_t n, RandTwister::ResultType *a) {
        const auto leftover(n & 0x7ULL);
        n >>= 3;
        switch(leftover) {
            case 0: do {
                    *a++ = twister_();
            case 7: *a++ = twister_();
            case 6: *a++ = twister_();
            case 5: *a++ = twister_();
            case 4: *a++ = twister_();
            case 3: *a++ = twister_();
            case 2: *a++ = twister_();
            case 1: *a++ = twister_();
                       } while(n--);
        }
    }
    void reseed(ResultType seed) {twister_.seed(seed);}
};

struct ThreadsafeRandTwister: public RandTwister {
    tthread::fast_mutex lock_;
    auto operator()() {
        lock_.lock();
        const auto ret(twister_());
        lock_.unlock();
        return ret;
    }
    // Generate a large number of random integers.
    auto operator()(size_t n, ResultType *a) {
        lock_.lock();
        RandTwister::operator ()(n, a);
        lock_.unlock();
    }
};

static RandTwister random_twist(std::rand());

template<typename T>
T randf() {return random_twist() * RandTwister::MAX_INV;}

} // namespace rand

#endif // #ifndef _EMP_RAND_H__
