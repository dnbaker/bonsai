#ifndef BONSAI_POLYHASH_H__
#define BONSAI_POLYHASH_H__
#include <array>
#include <cstdint>
#include <type_traits>

namespace bns {

/*
 * Rolling polynomial (Rabin-Karp) hash over GF(p), p = 2^61 - 1, with the
 * same interface as CyclicHash (eat, update, reverse_update, reset, seed and
 * a public hashvalue).
 *
 * For an n-gram c_0 .. c_{n-1}, each lane computes
 *     h = sum_i T[c_i] * B^(n - 1 - i) mod p
 * with a random base B and random character keys T. Unlike a cyclic hash in a
 * w-bit word, where positions i and i + w get the same rotation, distinct
 * positions always get distinct powers of B, so the hash has no period in n.
 * Two distinct n-grams collide in a lane with probability at most (n - 1) / p
 * over the choice of B.
 *
 * hashvalue is a bijective 64-bit mix of each lane; 64-bit output uses one
 * lane and 128-bit output two independent lanes.
 */
template<typename hashvaluetype=uint64_t, typename chartype=unsigned char>
class PolyHash {
    static constexpr uint64_t P = (uint64_t(1) << 61) - 1;
    static constexpr int NLANES = sizeof(hashvaluetype) > 8 ? 2: 1;
    static_assert(sizeof(hashvaluetype) == 8 || sizeof(hashvaluetype) == 16, "PolyHash produces 64- or 128-bit values");
    struct Lane {
        uint64_t base, base_inv, top; // B, B^-1 and B^(n-1) mod p
        std::array<uint64_t, 256> keys;
        uint64_t state;
    };
    int n_;
    std::array<Lane, NLANES> lanes_;

    static uint64_t reduce(unsigned __int128 x) {
        uint64_t r = uint64_t(x & P) + uint64_t(x >> 61);
        r = (r & P) + (r >> 61);
        return r >= P ? r - P: r;
    }
    static uint64_t mulmod(uint64_t a, uint64_t b) {return reduce((unsigned __int128)a * b);}
    static uint64_t addmod(uint64_t a, uint64_t b) {a += b; return a >= P ? a - P: a;}
    static uint64_t submod(uint64_t a, uint64_t b) {return a >= b ? a - b: a + P - b;}
    static uint64_t powmod(uint64_t a, uint64_t e) {
        uint64_t r = 1;
        for(; e; e >>= 1, a = mulmod(a, a))
            if(e & 1) r = mulmod(r, a);
        return r;
    }
    static uint64_t splitmix(uint64_t &x) {
        uint64_t z = (x += 0x9e3779b97f4a7c15ULL);
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        return z ^ (z >> 31);
    }
    // Uniform-enough draw from [lo, P).
    static uint64_t draw(uint64_t &x, uint64_t lo) {
        uint64_t v;
        do v = splitmix(x) >> 3; while(v >= P || v < lo);
        return v;
    }
    static uint64_t mix(uint64_t z) {
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        return z ^ (z >> 31);
    }
    void set_output() {
        if constexpr(NLANES == 1) {
            hashvalue = hashvaluetype(mix(lanes_[0].state));
        } else {
            hashvalue = (hashvaluetype(mix(lanes_[0].state)) << 64) | hashvaluetype(mix(lanes_[1].state ^ 0x5851f42d4c957f2dULL));
        }
    }
public:
    hashvaluetype hashvalue;

    // n is the n-gram length; the second argument (the word size of
    // CyclicHash) is accepted for interface compatibility and ignored.
    PolyHash(int n, int=0): n_(n) {seed(0);}

    void seed(uint64_t s1, uint64_t s2=0) {
        uint64_t x = s1 ^ (s2 * 0x9e3779b97f4a7c15ULL);
        for(auto &l: lanes_) {
            l.base = draw(x, uint64_t(1) << 32);
            l.base_inv = powmod(l.base, P - 2);
            l.top = powmod(l.base, n_ > 0 ? n_ - 1: 0);
            for(auto &k: l.keys) k = draw(x, 1);
        }
        reset();
    }
    void reset() {
        for(auto &l: lanes_) l.state = 0;
        set_output();
    }
    // Appends inchar: hash of ABC becomes hash of ABC[inchar].
    void eat(chartype inchar) {
        for(auto &l: lanes_)
            l.state = addmod(mulmod(l.state, l.base), l.keys[static_cast<uint8_t>(inchar)]);
        set_output();
    }
    // Hash of [outchar]ABC becomes hash of ABC[inchar].
    void update(chartype outchar, chartype inchar) {
        for(auto &l: lanes_) {
            const uint64_t rest = submod(l.state, mulmod(l.keys[static_cast<uint8_t>(outchar)], l.top));
            l.state = addmod(mulmod(rest, l.base), l.keys[static_cast<uint8_t>(inchar)]);
        }
        set_output();
    }
    // Hash of ABC[inchar] becomes hash of [outchar]ABC.
    void reverse_update(chartype outchar, chartype inchar) {
        for(auto &l: lanes_) {
            const uint64_t rest = mulmod(submod(l.state, l.keys[static_cast<uint8_t>(inchar)]), l.base_inv);
            l.state = addmod(rest, mulmod(l.keys[static_cast<uint8_t>(outchar)], l.top));
        }
        set_output();
    }
};

} // namespace bns

#endif /* BONSAI_POLYHASH_H__ */
