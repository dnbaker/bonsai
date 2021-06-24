#ifndef _EMP_SPACE_H__
#define _EMP_SPACE_H__

#include <vector>
#include <string>
#include <algorithm>
#include "kmerutil.h"

namespace bns {
using std::uint16_t;

using spvec_t = std::vector<uint16_t>;

inline u32 comb_size(const spvec_t &spaces) {
    u32 ret(spaces.size() + 1); // Since there's 1 fewer entry in spaces
    // We then increment the size of our comb for each space.
    for(const auto i: spaces) ret += i;
    return ret;
}

inline ks::string str(const spvec_t &vec) {
    ks::string ret;
    ret.resize(vec.size() << 1);
    for(const auto &val: vec) ret.sprintf("%u,", val);
    ret.pop();
    return ret;
}

static spvec_t parse_spacing(const char *ss, unsigned k) {
    if(!ss || *ss == '\0') return spvec_t(k - 1, 0); // No spaces

    spvec_t ret;
    while(ss && *ss) {
        int j(atoi(ss));
        ret.emplace_back(j);

        if(strchr(ss, 'x')) {
            ss = strchr(ss, 'x') + 1;
            for(int k(atoi(ss) - 1); k; k--) ret.emplace_back(j);
        }
        ss = strchr(ss, ',') + 1;
    }
    return ret;
}

struct Spacer {
    //static constexpr u32 max_k = sizeof(uint64_t) * CHAR_BIT / 2;

    // Instance variables
    spvec_t   s_; // Spaces to skip
    u32 k_; // Kmer size
    u32 c_; // comb size
    u32 w_; // window size

public:
    Spacer(unsigned k, uint32_t w, spvec_t spaces=spvec_t{}):
      s_(spaces.size() ? spaces: spvec_t(k - 1, 0)),
      k_(k),
      c_(comb_size(s_)),
      w_(std::max((int)c_, (int)w))
    {
        //if(k > max_k) LOG_WARNING("Provided k %u greater than can uniquely be described by 64-bit integers (%u).\n", k_, max_k);
        for(auto &i: s_) ++i; // Convert differences into offsets
        if(s_.size() + 1 != k) {
            LOG_EXIT("Error: input vector must have size 1 less than k. k: %u. size: %zu.\n",
                     k, s_.size());
        }
    }
    u32 k() const {return k_;}
    u32 w() const {return w_;}
    u32 c() const {return c_;}
    spvec_t spaces() const {return s_;}
    Spacer& operator=(const Spacer &o) = default;
    void resize(unsigned k, unsigned w, std::string space_string="") {
        *this = Spacer(k, w, space_string.data());
    }
    Spacer(unsigned k, uint32_t w, const char *space_string):
        Spacer(k, w, parse_spacing(space_string, k)) {}
    bool unspaced() const {
        return std::find_if(s_.begin(), s_.end(), [](auto x) {return x != 1;}) == s_.end();
    }
    bool unwindowed() const {
        return k_ == w_;
    }
    Spacer(unsigned k): Spacer(k, k) {}
    Spacer(const Spacer &other): s_(other.s_), k_(other.k_), c_(other.c_), w_(other.w_) {}
    auto write(u128 kmer, std::FILE *fp=stdout) const {
        char static_buf[256];
        char *buf = c_ <= sizeof(static_buf) ? static_buf: static_cast<char *>(std::malloc(c_));
        char *bp = buf;
        auto offset = ((k_ - 1) * 2);
        auto it = s_.begin();
        *bp++ = num2nuc((kmer >> offset) & 0x3u);
        do {
            offset -= 2;
            for(int i = *it++; i-- > 1; *bp++ = '-');
            *bp++ = num2nuc((kmer >> offset) & 0x3u);
        } while(it != s_.end());

        const int ret = std::fwrite(buf, 1, c_, fp);
        if(c_ > sizeof(static_buf)) std::free(buf);
        return ret;
    }
    auto write(u64 kmer, std::FILE *fp=stdout) const {
        char static_buf[256];
        char *buf = c_ <= sizeof(static_buf) ? static_buf: static_cast<char *>(std::malloc(c_));
        char *bp = buf;
        auto offset = ((k_ - 1) * 2);
        auto it = s_.begin();
        *bp++ = num2nuc((kmer >> offset) & 0x3u);
        do {
            offset -= 2;
            for(int i = *it++; i-- > 1; *bp++ = '-');
            *bp++ = num2nuc((kmer >> offset) & 0x3u);
        } while(it != s_.end());

        const int ret = std::fwrite(buf, 1, c_, fp);
        if(c_ > sizeof(static_buf)) std::free(buf);
        return ret;
    }
    std::string to_string(u128 kmer) const {
        std::string ret;
        ret.reserve(c_ - k_ + 1);
        int offset = ((k_ - 1) << 1);
        ret.push_back(num2nuc((kmer >> offset) & 0x3u));
        for(auto s: s_) {
            assert(offset >= 0);
            offset -= 2;
            while(s-->1) ret.push_back('-');
            ret.push_back(num2nuc((kmer >> offset) & 0x3u));
        }
        return ret;
    }
    std::string to_string(u64 kmer) const {
        std::string ret;
        ret.reserve(c_ - k_ + 1);
        int offset = ((k_ - 1) << 1);
        ret.push_back(num2nuc((kmer >> offset) & 0x3u));
        for(auto s: s_) {
            assert(offset >= 0);
            offset -= 2;
            while(s-->1) ret.push_back('-');
            ret.push_back(num2nuc((kmer >> offset) & 0x3u));
        }
        return ret;
    }
    ~Spacer() {}
    spvec_t sub1() const {spvec_t ret(s_);std::transform(ret.begin(), ret.end(), ret.begin(), [](auto x) {return x - 1;}); return ret;}
};

} // namespace bns

#endif // #ifndef _EMP_SPACE_H__
