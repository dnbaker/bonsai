#pragma once
#include <type_traits>
#include <stdexcept>
#include <cstdlib>

namespace kmp {
using std::size_t;

template<typename T>
struct DestructIf {
    void operator()(const T *v) const {const_cast<T *>(v)->~T();}
};

template<typename T, typename FreeFunc=DestructIf<T>>
class Pool {
    // Pointer memory pool
    size_t cnt_, n_, max_;
    T **buf_;
    FreeFunc func_;
public:
    Pool(): cnt_(0), n_(0), max_(0), buf_(0), func_() {}
    T *malloc() {
        ++cnt_;
        if(n_ == 0) return static_cast<T *>(std::malloc(sizeof(T)));
        return buf_[--n_];
    }
    T *calloc() {
        ++cnt_;
        if(n_ == 0) return static_cast<T *>(std::calloc(1, sizeof(T)));
        auto ret = buf_[--n_];
        std::memset(ret, 0, sizeof(T)); // Zero memory pointed to.
        return ret;
    }
    template<typename... Args>
    T *placement_new(Args &&... args) {
        ++cnt_;
        if(n_ == 0) return new T(std::forward<Args>(args)...);
        T *ret = buf_[n_];
        func_(ret); // This should be its destructor.
        return buf_[--n_];
    }
    void free(T *p) {
        --cnt_;
        if(n_ == max_) {
            max_ = max_ ? max_ << 1: 16;
            if((buf_ = static_cast<T **>(std::realloc(buf_, sizeof(T *) * max_))) == nullptr) throw std::bad_alloc();
        }
        buf_[n_++] = p;
    }
    ~Pool() {
        for(size_t k = 0; k < n_; ++k) {
            func_(buf_[k]);
            std::free(buf_[k]);
        }
        std::free(buf_);
    }
};

} // namespace kmp
