#pragma once
#include <new>         // For placement new
#include <cassert>     // For assert
#include <cstdlib>     // For std::size_t
#include <cstdint>     // For std::uint*_ts
#include <stdexcept>   // For std::bad_alloc
#include <string>      // for std::string [for exception handling]
#include <type_traits> // For std::enable_if_t/std::is_unsigned_v

namespace circ {
using std::size_t;
using std::uint64_t;
using std::uint32_t;
using std::uint16_t;
using std::uint8_t;

static inline uint64_t roundup(uint64_t x) {
    --x, x|=x>>1, x|=x>>2, x|=x>>4, x|=x>>8, x|=x>>16, x|=x>>32;
    return ++x;
}
static inline uint32_t roundup(uint32_t x) {
    --x, x|=x>>1, x|=x>>2, x|=x>>4, x|=x>>8, x|=x>>16;
    return ++x;
}
static inline uint16_t roundup(uint16_t x) {
    --x, x|=x>>1, x|=x>>2, x|=x>>4, x|=x>>8;
    return ++x;
}
static inline uint8_t roundup(uint8_t x) {
    --x, x|=x>>1, x|=x>>2, x|=x>>4;
    return ++x;
}

template<typename T, typename SizeType=uint32_t, typename=std::enable_if_t<std::is_unsigned_v<T>>>
class FastCircularQueue {
    // A circular queue in which extra memory has been allocated up to a power of two.
    // This allows us to use bitmasks instead of modulus operations.
    // This circular queue is NOT threadsafe. Its purpose is creating a double-ended queue without
    // the overhead of a doubly-linked list.
    const SizeType mask_;
    SizeType      start_;
    SizeType       stop_;
    T             *data_;

public:
    using size_type = SizeType;
    class iterator {
        // TODO: increment by an integral quantity.
        FastCircularQueue &ref_;
        SizeType           pos_;
    public:
        iterator(FastCircularQueue &ref, SizeType pos): ref_(ref), pos_(pos) {}
        iterator(const iterator &other): ref_(other.ref_), pos_(other.pos_) {}
        T &operator*() {
            return ref_.data_[pos_];
        }
        const T &operator*() const noexcept {
            return ref_.data_[pos_];
        }
        iterator &operator++() noexcept {
            ++pos_;
            pos_ &= ref_.mask_;
            return *this;
        }
        iterator operator++(int) noexcept {
            iterator copy(*this);
            this->operator++();
            return copy;
        }
        bool operator==(const iterator &other) const noexcept {
            return pos_ == other.pos_;
        }
        bool operator!=(const iterator &other) const noexcept {
            return pos_ != other.pos_;
        }
        bool operator<(const iterator &other) const noexcept {
            return pos_ < other.pos_;
        }
        bool operator<=(const iterator &other) const noexcept {
            return pos_ <= other.pos_;
        }
        bool operator<(const iterator &other) const noexcept {
            return pos_ < other.pos_;
        }
        bool operator>=(const iterator &other) const noexcept {
            return pos_ >= other.pos_;
        }
    };
    class const_iterator {
        const FastCircularQueue &ref_;
        SizeType                 pos_;
    public:
        const_iterator(const FastCircularQueue &ref, SizeType pos): ref_(ref), pos_(pos) {}
        const_iterator(const const_iterator &other): ref_(other.ref_), pos_(other.pos_) {}
        const T &operator*() const noexcept {
            return ref_.data_[pos_];
        }
        const_iterator &operator++() noexcept {
            ++pos_;
            pos_ &= ref_.mask_;
            return *this;
        }
        const_iterator operator++(int) noexcept {
            const_iterator copy(*this);
            this->operator++();
            return copy;
        }
        bool operator==(const const_iterator &other) const noexcept {
            return pos_ == other.pos_;
        }
        bool operator!=(const const_iterator &other) const noexcept {
            return pos_ != other.pos_;
        }
        bool operator<(const const_iterator &other) const noexcept {
            return pos_ < other.pos_;
        }
        bool operator<=(const const_iterator &other) const noexcept {
            return pos_ <= other.pos_;
        }
        bool operator>(const const_iterator &other) const noexcept {
            return pos_ > other.pos_;
        }
        bool operator>=(const const_iterator &other) const noexcept {
            return pos_ >= other.pos_;
        }
    };
    iterator begin() noexcept {
        return iterator(*this, start_);
    }
    iterator end() noexcept {
        return iterator(*this, stop_);
    }
    const_iterator cbegin() const noexcept {
        return const_iterator(*this, start_);
    }
    const_iterator cend()   const noexcept {
        return const_iterator(*this, stop_);
    }
    FastCircularQueue(SizeType size):
            mask_(roundup(size) - 1),
            start_(0), stop_(0),
            data_(static_cast<T *>(std::malloc(sizeof(T) * capacity())))
    {
        assert((mask_ & (mask_ - 1)) == 0);
        if(data_ == nullptr) throw std::bad_alloc();
    }
    template<typename... Args>
    T &push(Args &&... args) {
        size_type ind = stop_;
        ++stop_; stop_ &= mask_;
        if(__builtin_expect(stop_ == start_, 0)) throw std::runtime_error(std::string("Overfull buffer of size ") + std::to_string(capacity()) + ". Abort!");
        return *(new(data_ + ind) T(std::forward<Args>(args)...));
    }
    T pop() {
        if(__builtin_expect(stop_ == start_, 0)) throw std::runtime_error("Popping item from empty buffer. Abort!");
        T ret(std::move(data_[start_++]));
        start_ &= mask_;
        return ret; // If unused, the std::move causes it to leave scope and therefore be destroyed.
    }
    template<typename... Args>
    T push_pop(Args &&... args) {
        T ret(pop());
        push(std::forward<Args>(args)...);
        return ret;
    }
    ~FastCircularQueue() noexcept {
        for(size_type i(start_); i != stop_; data_[i++].~T(), i &= mask_);
        std::free(data_);
    }
    size_type capacity() const noexcept {return mask_ + 1;}
    size_type size()     const noexcept {return (stop_ - start_) & mask_;}
    
} // FastCircularQueue

} // namespace circ
