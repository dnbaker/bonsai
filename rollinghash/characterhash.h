#ifndef CHARACTERHASH
#define CHARACTERHASH

typedef unsigned long long uint64;
typedef unsigned int uint32;
typedef unsigned int uint;

#include <cassert>
#include <iostream>
#include <stdexcept>
#include <random>
using std::cerr;
using std::runtime_error;
using std::endl;

template <typename hashvaluetype>
#if __cplusplus >= 201402L
constexpr
#endif
hashvaluetype maskfnc(int bits) {
    assert(bits>0);
    assert(unsigned(bits)<=sizeof(hashvaluetype)*8);
    hashvaluetype x = static_cast<hashvaluetype>(1) << (bits - 1);
    return x ^ (x - 1);
}
template<typename T>
inline T roundup(T v) {
    if(sizeof(T) == 4) return roundup(uint32_t(v));
    if(sizeof(T) == 8) return roundup(uint64_t(v));
    if(sizeof(T) == 16) return roundup(__uint128_t(v));
    return 0;
}
template<> inline uint32_t roundup(uint32_t v) {
    --v;
    v|= v>>1;
    v|= v>>2;
    v|= v>>4;
    v|= v>>8;
    v|= v>>16;
    return ++v;
}
template<> inline uint64_t roundup(uint64_t v) {
    --v;
    v|= v>>1;
    v|= v>>2;
    v|= v>>4;
    v|= v>>8;
    v|= v>>16;
    v|= v>>32;
    return ++v;
}
template<> inline __uint128_t roundup(__uint128_t v) {
    --v;
    v|= v>>1;
    v|= v>>2;
    v|= v>>4;
    v|= v>>8;
    v|= v>>16;
    v|= v>>32;
    v|= v>>64;
    return ++v;
}
template<> inline int32_t roundup(int32_t v) {return roundup(uint32_t(v));}
template<> inline int64_t roundup(int64_t v) {return roundup(uint64_t(v));}
template<> inline __int128_t roundup(__int128_t v) {return roundup(__uint128_t(v));}

template<typename T>
inline void hack(T &x, std::mt19937_64 &rng) {x = rng();}
template<> inline void hack(__uint128_t &x, std::mt19937_64 &rng) {
    x = rng(); x <<= 64; x |= rng();
}

template <typename hashvaluetype = uint32, typename chartype =  unsigned char>
class CharacterHash {
public:

    static constexpr size_t HVSIZE = sizeof(hashvaluetype);
    static_assert(HVSIZE == 2 || HVSIZE == 1 || HVSIZE == 4 || HVSIZE == 8 || HVSIZE == 16, "Must be a supported type");
    CharacterHash(hashvaluetype maxval, uint32 seed1=0, uint32 seed2=0x1337) {
        if(seed1 == 0) seed1 = std::rand();
        seed(maxval, seed1, seed2);
    }
    void seed(hashvaluetype maxval, uint32_t seed1, uint64_t seed2=0) {
        std::mt19937_64 randomgenerator(seed1), randomgeneratorbase(seed2);
        randomgenerator.seed(seed1);
        hashvaluetype tmaxval = roundup(maxval) - 1;
        for(size_t k =0; k<nbrofchars; ++k) {
            hashvaluetype next;
            do { hack(next, randomgenerator); next &= tmaxval;} while(next > maxval);
            hashvalues[k] = next;
        }
    }

    enum {nbrofchars = 1 << ( sizeof(chartype)*8 )};

    hashvaluetype hashvalues[nbrofchars];
};

#endif
