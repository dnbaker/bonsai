#ifndef _LINEAR_SET_H__
#define _LINEAR_SET_H__
#include <type_traits>
#include <stdexcept>
#include <cstring>
#include <ratio>
#include <memory>
#include <cstdlib>

namespace del {
struct free_deleter{
    template <typename T>
    void operator()(T *p) const {
        std::free(const_cast<std::remove_const_t<T>*>(p));
    }
};
}

namespace linear {

template<typename T, typename SizeType=size_t, class ResizeRatio=std::ratio<3,2>>
class set {
    SizeType n_, m_;
    std::unique_ptr<T, del::free_deleter> data_;

    static constexpr const double RESIZE_FACTOR = static_cast<double>(ResizeRatio::num) / static_cast<double>(ResizeRatio::den);

public:
    using size_type  = SizeType;
    using value_type = T;
    using difference_type = std::make_signed_t<size_type>;

    template<typename... FactoryArgs>
    set(size_type n=0, bool init=!std::is_pod<T>::value, FactoryArgs &&...args): n_{n}, m_{n}, data_{static_cast<T *>(std::malloc(sizeof(T) * n))} {
        if(data_ == nullptr)
            throw std::bad_alloc();
        if (init)
            for(size_type i(0); i < n_; ++i)
                new(data_.get() + i) T(std::forward<FactoryArgs>(args)...);
    }
    set(set &&other): n_(0), m_(0), data_(nullptr) {
        std::swap(*this, other);
    }
    set(const set &other): n_(other.n_), m_(other.m_), data_(std::malloc(sizeof(T) * n_)) {
        if(data_ == nullptr)
            throw std::bad_alloc();
        std::copy(other.begin(), other.end(), begin());
    }
    const T *cbegin() const {return data_;}
    const T *cend()   const {return data_.get() + n_;}
    const T *begin()  const {return data_;}
    const T *end()    const {return data_.get() + n_;}
    T *begin()        {return data_.get();}
    T *end()          {return data_.get() + n_;}
    T *data()         {return data_.get();}
    const T *data() const {return data_;}
    T &back() {return data_.get()[n_ - 1];}
    const T &back() const {return data_.get()[n_ - 1];}
    T &front() {return data_.get()[0];}
    const T &front() const {return data_.get()[0];}
    auto size() const {return n_;}
    auto capacity() const {return m_;}
    auto find(const T &val) {
        return std::find(begin(), end(), val);
    }
    template<typename Predicate>
    auto find_if(const T &val, Predicate p) {
        return std::find_if(begin(), end(), p);
    }
    bool contains(const T &val) {
        return find(val) != end();
    }
    T *insert(const T &val) {
        auto it(find(val));
        if(it == end()) {
            it = &push_back(val);
        }
        return it;
    }
    // Note: push_back and emplace_back are the same *except* that push_back changes size multiplicatively.
    template<typename... Args>
    auto &emplace_back(Args&& ... args) {
        if(n_ + 1 > m_) {
            m_ = std::max(static_cast<size_type>(m_ * RESIZE_FACTOR), m_ + 1);
            auto tmp(static_cast<T *>(std::realloc(data_, sizeof(T) * (m_))));
            if(tmp == nullptr) throw  std::bad_alloc();
            data_.release();
            data_.reset(tmp);
        }
        new(data_ + n_++) T(std::forward<Args>(args)...);
        return back();
    }
    template<typename... Args>
    auto &emplace_front(Args&& ... args) {
        if(n_ + 1 > m_) {
            m_ = std::max(static_cast<size_type>(m_ * RESIZE_FACTOR), m_ + 1);
            auto tmp(static_cast<T *>(std::realloc(data_, sizeof(T) * (m_))));
            if(tmp == nullptr) throw  std::bad_alloc();
            data_.release();
            data_.reset(tmp);
        }
        std::memmove(data_.get() + sizeof(T), data_, n_ * sizeof(T));
        new(data_) T(std::forward<Args>(args)...);
        ++n_;
        return front();
    }
    auto &push_back(const T &val) {
        if(n_ + 1 > m_) {
            m_ = std::max(static_cast<size_type>(m_ * RESIZE_FACTOR), m_ + 1);
            auto tmp(static_cast<T *>(std::realloc(data_.get(), sizeof(T) * (m_))));
            if(tmp == nullptr) throw std::bad_alloc();
            data_.release();
            data_.reset(tmp);
        }
        data_.get()[n_++] = val;
        return back();
    }
    
    void zero() {std::memset(data_.get(), 0, sizeof(T) * n_);} // DOES NOT CALL DESTRUCTORS
    void reserve(size_t newsize) {
        if(newsize > m_) {
            auto tmp(static_cast<T*>(std::realloc(data_.get(), sizeof(T) * newsize)));
            if(tmp == nullptr) throw std::bad_alloc();
            data_.release();
            data_.reset(tmp);
        }
    }
    template<typename... Args>
    void resize(size_t newsize, bool init=true, Args &&...args) {
        reserve(newsize);
        if(init) {
            for(;n_ < m_;) new(data_.get() + n_++) T(std::forward<Args>(args)...);
        }
    }
    void shrink_to_fit() {
        if(m_ > n_) {
            auto tmp(static_cast<T*>(std::realloc(data_, sizeof(T) * n_)));
            if(tmp == nullptr) throw std::bad_alloc();
            data_.release();
            data_.reset(tmp);
            m_ = n_;
        }
    }
    T &operator[](size_type idx) {return data_[idx];}
    const T &operator[](size_type idx) const {return data_[idx];}
    ~set(){
        if constexpr(!std::is_trivially_destructible<T>::value) {
            for(auto &el: *this) el.~T();
        }
    }
};

}

#endif // #ifndef _LINEAR_SET_H__
