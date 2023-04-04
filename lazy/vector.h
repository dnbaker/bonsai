#ifndef _LZY_VEC_H__
#define _LZY_VEC_H__
#include <type_traits>
#include <stdexcept>
#include <cstring>
#include <cstdlib>
#include <memory>
#include <algorithm>

#ifndef LAZY_PUSH_BACK_RESIZING_FACTOR
#define LAZY_PUSH_BACK_RESIZING_FACTOR 1.25
#endif


namespace lazy {

enum initialization: bool {
    LAZY_VEC_NOINIT = false,
    LAZY_VEC_INIT = true
};

template<typename T, typename size_type=uint32_t>
class vector {
    size_type n_, m_;
    T *data_;
    static constexpr double PUSH_BACK_RESIZING_FACTOR = LAZY_PUSH_BACK_RESIZING_FACTOR;

public:
    using value_type = T;
    using const_pointer_type = const T *;
    using pointer_type = T *;
    template<typename W>
    struct is_pod: public std::integral_constant<bool, std::is_trivial<W>::value && std::is_standard_layout<W>::value> {};
    template<typename... FactoryArgs>
    vector(size_type n=0, bool init=!is_pod<T>::value, FactoryArgs &&...args): n_{n}, m_{n}, data_{n ? static_cast<pointer_type>(std::malloc(sizeof(T) * n)): nullptr} {
        if (init)
            for(size_type i(0); i < n_; ++i)
                new(data_ + i) T(std::forward<FactoryArgs>(args)...);
    }
    vector(const vector &other): n_(other.n_), m_(other.n_), data_{static_cast<pointer_type>(std::malloc(sizeof(T) * m_))} {
        std::copy(std::cbegin(other), std::cend(other), begin());
    }
    vector(vector &&other): n_(other.n_), m_(other.m_), data_(other.data_) {
        other.m_ = other.n_ = 0;
        other.data_ = nullptr;
    }
    vector(std::initializer_list<T> il): n_(il.size()), m_(il.size()), data_{n_ ? static_cast<pointer_type>(std::malloc(sizeof(T) * n_)): nullptr} {
        std::move(il.begin(), il.end(), std::begin(*this));
    }
    const_pointer_type cbegin() const {return data_;}
    const_pointer_type cend()   const {return data_ + n_;}
    const_pointer_type begin()  const {return data_;}
    const_pointer_type end()    const {return data_ + n_;}
    const_pointer_type data()   const {return data_;}

    pointer_type data()         {return data_;}

    pointer_type begin()        {return data_;}
    pointer_type end()          {return data_ + n_;}
    auto &back() {return data_[n_ - 1];}
    const auto &back() const {return data_[n_ - 1];}
    auto &front() {return *data_;}
    const auto &front() const {return *data_;}
    auto size() const {return n_;}
    auto capacity() const {return m_;}
    vector &operator=(const vector &o) {
        if(m_ < o.m_) data_ = static_cast<pointer_type>(std::realloc(data_, sizeof(T) * o.m_));
        n_ = o.n_;
        m_ = o.m_;
        std::copy(o.begin(), o.end(), begin());
        return *this;
    }
    vector &operator=(vector &&o) {
        n_ = o.n_;
        m_ = o.m_;
        data_ = o.data_;
        o.release();
        return *this;
    }
    template<typename osize_type>
    bool operator==(const ::lazy::vector<T, osize_type> &other) const {
        if(size() != other.size()) return false;
        for(typename std::common_type_t<osize_type, size_type>i(0);i < n_; ++i)
            if(data_[i] != other[i])
                return false;
        return true;
    }
    // Note: push_back and bnslace_back are the same *except* that push_back changes size multiplicatively.
    template<typename... Args>
    auto &emplace_back(Args&& ... args) {
        if(n_ + 1 > m_) {
            data_ = static_cast<pointer_type>(std::realloc(data_, sizeof(T) * (n_ + 1)));
#if !NDEBUG
            if(__builtin_expect(m_ + 1 < m_, 0)) {
                char buf[256];
                std::sprintf(buf, "Type of %zu bytes is not enough to hold the full values. (%zu, %zu)\n",
                            sizeof(m_), m_, m_ + 1);
                throw std::runtime_error(buf);
            }
#else
            ++m_;
#endif
        }
        new(data_ + n_++) T(std::forward<Args>(args)...);
        return back();
    }
    template<typename... Args>
    auto &push_back(Args&& ... args) {
        if(n_ + 1 > m_)
            reserve(std::max(static_cast<size_type>(m_ * PUSH_BACK_RESIZING_FACTOR),
                             m_ + 1));
        new(data_ + n_++) T(std::forward<Args>(args)...);
        return back();
    }
    void zero() {std::memset(data_, 0, sizeof(T) * n_);}
    void reserve(size_t newsize) {
        if(newsize > m_) {
            auto tmp(static_cast<T*>(std::realloc(data_, sizeof(T) * newsize)));
            if(tmp == nullptr) throw std::bad_alloc();
            data_ = tmp;
        }
    }
    template<typename... Args>
    void resize(size_t newsize, bool init=true, Args &&...args) {
        reserve(newsize);
        if(init) {
            for(;n_ < m_;) new(data_ + n_++) T(std::forward<Args>(args)...);
        }
    }
    void shrink_to_fit() {
        if(m_ > n_) {
            auto tmp(static_cast<T*>(std::realloc(data_, sizeof(T) * n_)));
            if(tmp == nullptr) throw std::bad_alloc();
            data_ = tmp;
            m_ = n_;
        }
    }
    T &operator[](size_type idx) {return data_[idx];}
    const T &operator[](size_type idx) const {return data_[idx];}
    ~vector(){
#if __cplusplus < 201703L
        if(!std::is_trivially_destructible<T>::value) {
#else
        if constexpr(!std::is_trivially_destructible_v<T>) {
#endif
            for(auto &el: *this) el.~T();
        }
        std::free(data_);
    }
    void release() {n_ = m_ = 0; data_ = 0;}
};

}

#endif // #ifndef _LZY_VEC_H__
