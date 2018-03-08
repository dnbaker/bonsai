#pragma once
#include <new>         // For placement new
#include <cassert>     // For assert
#include <cstdlib>     // For std::size_t
#include <cstdint>     // For std::uint*_ts
#include <stdexcept>   // For std::bad_alloc
#include <string>      // for std::string [for exception handling]
#include <type_traits> // For std::enable_if_t/std::is_unsigned_v
#include <cstring>     // For std::memcpy
#include <vector>      // For converting to vector

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

template<typename T, typename SizeType=uint32_t>
class FastCircularQueue {
    // A circular queue in which extra memory has been allocated up to a power of two.
    // This allows us to use bitmasks instead of modulus operations.
    // This circular queue is NOT threadsafe. Its purpose is creating a double-ended queue without
    // the overhead of a doubly-linked list.
    SizeType  mask_;
    SizeType start_;
    SizeType  stop_;
    T        *data_;
    static_assert(std::is_unsigned_v<SizeType>, "Must be unsigned");

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
        bool operator>(const iterator &other) const noexcept {
            return pos_ > other.pos_;
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
            mask_(roundup(size + 1) - 1),
            start_(0), stop_(0),
            data_(static_cast<T *>(std::malloc((mask_ + 1) * sizeof(T))))
    {
        assert((mask_ & (mask_ + 1)) == 0);
        if(data_ == nullptr) {
            throw std::bad_alloc();
        }
    }
    void resize(size_type new_size) {
        if(__builtin_expect(new_size < mask_, 0)) throw std::runtime_error("Attempting to resize to value smaller than queue's size, either from user error or overflowing the size_type. Abort!");
        new_size = roundup(new_size); // Is this necessary? We can hide resize from the user and then cut out this call.
        auto tmp = std::realloc(data_, new_size * sizeof(T));
        if(tmp == nullptr) throw std::bad_alloc();
        data_ = static_cast<T *>(tmp);
        if(start_ == stop_) {
            if(start_) {
                stop_ = mask_ + 1;
                auto tmp = static_cast<T *>(std::malloc((stop_) * sizeof(T)));
                if(tmp == nullptr) throw std::bad_alloc();
                std::memcpy(tmp, data_ + start_, (stop_ - start_) * sizeof(T));
                std::memcpy(tmp + (stop_ - start_), data_, stop_ * sizeof(T));
                std::memcpy(data_, tmp, (stop_) * sizeof(T));
                std::free(tmp);
                start_ = 0;
            }
        } else if(stop_ < start_) {
            auto tmp = static_cast<T *>(std::malloc((mask_ + 1) * sizeof(T)));
            if(tmp == nullptr) throw std::bad_alloc();
            std::memcpy(tmp, data_ + start_, ((mask_ + 1) - start_) * sizeof(T));
            std::memcpy(tmp + ((mask_ + 1) - start_), data_, stop_ * sizeof(T));
            auto sz = (stop_ - start_) & mask_;
            std::memcpy(data_, tmp, sz);
            std::free(tmp);
            start_ = 0;
            stop_ = sz;
        }
        mask_ = ((mask_ + 1) << 1) - 1;
    }
    // Does not yet implement push_front.
    template<typename... Args>
    T &push_back(Args &&... args) {
        if(__builtin_expect(((stop_ + 1) & mask_) == start_, 0)) resize((mask_ + 1) << 1);
        size_type ind = stop_;
        ++stop_; stop_ &= mask_;
        return *(new(data_ + ind) T(std::forward<Args>(args)...));
    }
    template<typename... Args>
    T &emplace_back(Args &&... args) {
        return push_back(std::forward<Args>(args)...); // Interface compatibility.
    }
    T pop() {
        if(__builtin_expect(stop_ == start_, 0)) throw std::runtime_error("Popping item from empty buffer. Abort!");
        T ret(std::move(data_[start_++]));
        start_ &= mask_;
        return ret; // If unused, the std::move causes it to leave scope and therefore be destroyed.
    }
    T pop_back() {
        if(__builtin_expect(stop_ == start_, 0)) throw std::runtime_error("Popping item from empty buffer. Abort!");
        T ret(std::move(data_[--stop_]));
        start_ &= mask_;
        return ret; // If unused, the std::move causes it to leave scope and therefore be destroyed.
    }
    T pop_front() {
        return pop(); // Interface compatibility with std::list.
    }
    template<typename... Args>
    T push_pop(Args &&... args) {
        T ret(pop());
        push(std::forward<Args>(args)...);
        return ret;
    }
    T &back() {
        return data_[stop_ - 1];
    }
    const T &back() const {
        return data_[stop_ - 1];
    }
    T &front() {
        return data_[start_];
    }
    const T &front() const {
        return data_[start_];
    }
    ~FastCircularQueue() {clear(); std::free(data_);}
    size_type capacity() const noexcept {return mask_;}
    size_type size()     const noexcept {return (stop_ - start_) & mask_;}
    std::vector<T> to_vector() const {
        std::vector<T> ret;
        ret.reserve(this->size());
        for(size_type i(start_); i != stop_; ret.emplace_back(std::move(data_[i])), ++i, i &= mask_);
        return ret;
    }
    void clear() {
        if constexpr(std::is_destructible_v<T>) {
            for(size_type i(start_); i != stop_; data_[i++].~T(), i &= mask_);
        }
        start_ = stop_ = 0;
    }
}; // FastCircularQueue
template<typename T, typename SizeType>
using deque = FastCircularQueue<T, SizeType>;

} // namespace circ
