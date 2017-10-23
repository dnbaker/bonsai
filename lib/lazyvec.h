#ifndef _LZY_VEC_H__
#define _LZY_VEC_H__
#include <type_traits>
#include <stdexcept>
#include <cstring>
#include <cstdlib>

namespace lazy {

template<typename T, typename size_type=std::size_t, bool init=!std::is_pod<T>::value>
class vector {
    T *data_;
    size_type n_, m_;
    static constexpr double PUSH_BACK_RESIZING_FACTOR = 1.25;

public:
    using value_type = T;
    template<typename... FactoryArgs>
    vector(size_t n=0, FactoryArgs &&...args): n_{n}, m_{n}, data_{static_cast<T *>(std::malloc(sizeof(T) * n))} {
        if constexpr(init) {
            for(size_type i(0); i < n_; ++i) new(data_ + i) T(std::forward<FactoryArgs>(args)...);
        }
    }
    const T *cbegin() const {return data_;}
    const T *cend()   const {return data_ + n_;}
    const T *begin()  const {return data_;}
    const T *end()          const {return data_ + n_;}
    T *begin()        {return data_;}
    T *end()          {return data_ + n_;}
    T *data()         {return data_;}
    const T *data() const {return data_;}
    auto &back() {return data_[n_ - 1];}
    const auto &back() const {return data_[n_ - 1];}
    auto &front() {return data_[0];}
    const auto &front() const {return data_[0];}
    auto size() const {return n_;}
    auto capacity() const {return m_;}
    template<typename osize_type, bool oinit>
    bool operator==(const ::lazy::vector<T, osize_type, oinit> &other) const {
        if(size() != other.size()) return false;
        using IterType = std::common_type_t<osize_type, size_type>;
        for(IterType i(0); i < n_; ++i) if(data_[i] != other[i]) return false;
        return true;
    }
    // Note: push_back and emplace_back are the same *except* that push_back changes size multiplicatively.
    template<typename... Args>
    auto &emplace_back(Args&& ... args) {
        if(n_ + 1 > m_) {
            data_ = static_cast<T *>(std::realloc(data_, sizeof(T) * (n_ + 1)));
            ++m_;
        }
        new(data_ + n_++) T(std::forward<Args>(args)...);
        return back();
    }
    template<typename... Args>
    auto &push_back(Args&& ... args) {
        if(n_ + 1 > m_) {
            m_ *= PUSH_BACK_RESIZING_FACTOR;
            data_ = static_cast<T *>(std::realloc(data_, sizeof(T) * (m_)));
        }
        new(data_ + n_++) T(std::forward<Args>(args)...);
        return back();
    }
    void zero() {std::memset(data_, 0, sizeof(T) * n_);}
    void reserve(std::size_t newsize) {
        if(newsize > m_) {
            auto tmp(static_cast<T*>(std::realloc(data_, sizeof(T) * newsize)));
            if(tmp == nullptr) throw std::bad_alloc();
            data_ = tmp;
        }
    }
    template<typename... Args>
    void resize(std::size_t newsize, bool initialize=true, Args &&...args) {
        reserve(newsize);
        if(initialize) {
            for(;n_ < m_;) new(data_ + n_++) T(std::forward<Args>(args)...);
        }
    }
    T &operator[](size_type idx) {return data_[idx];}
    const T &operator[](size_type idx) const {return data_[idx];}
    ~vector(){
        if constexpr(!std::is_trivially_destructible<T>::value) {
            for(auto &el: *this) el.~T();
        }
        std::free(data_);
    }
};

}

#endif // #ifndef _LZY_VEC_H__
