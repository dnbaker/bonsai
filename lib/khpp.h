#ifndef KHPP_H__
#define KHPP_H__
#include <functional>
#include <shared_mutex>
#include <mutex>
#include <cstring>
#include "tinythreadpp/source/fast_mutex.h"
using std::shared_mutex;

#if __GNUC__ >= 7
#  define CONSTEXPR_IF if constexpr
#else
#  define CONSTEXPR_IF if
#endif

#define PACKED

#ifndef PACKED
#  if __GNUC__ || __clang__
#    define PACKED __attribute__((packed))
#  else
#    define PACKED
#  endif
#endif

#ifndef kroundup64
#define kroundup64(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, (x)|=(x)>>32, ++(x))
#endif

namespace kh {
using u64 = std::uint64_t;

template<typename K, typename V>
struct packed_pair {
    K first;
    V second;
    packed_pair(K &&k, V &&v): first(k), second(v) {}
    bool operator==(const packed_pair &other) {
        return std::tie(first, second) == std::tie(other.first, other.second);
    }
    static constexpr bool trivially_destructible() {
        return std::is_trivially_destructible<K>::value && std::is_trivially_destructible<V>::value;
    }
} PACKED;

template<typename K, typename V, bool is_map>
struct khel_t;

template<typename K, typename V>
struct khel_t<K, V, true> {
    using Type = packed_pair<K, V>;
};

template<typename K, typename V>
struct khel_t<K, V, false> {
    using Type = K;
};

#define __ac_isempty(flag, i) ((flag[i>>4]>>((i&0xfU)<<1))&2)
#define __ac_isdel(flag, i) ((flag[i>>4]>>((i&0xfU)<<1))&1)
#define __ac_iseither(flag, i) ((flag[i>>4]>>((i&0xfU)<<1))&3)
#define __ac_set_isdel_false(flag, i) (flag[i>>4]&=~(1ull<<((i&0xfU)<<1)))
#define __ac_set_isempty_false(flag, i) (flag[i>>4]&=~(2ull<<((i&0xfU)<<1)))
#define __ac_set_isboth_false(flag, i) (flag[i>>4]&=~(3ull<<((i&0xfU)<<1)))
#define __ac_set_isdel_true(flag, i) (flag[i>>4]|=1ull<<((i&0xfU)<<1))

template<typename khkey_t, typename khval_t,
         typename hash_func=std::hash<khkey_t>, class hash_equal=std::equal_to<khkey_t>,
         bool is_map=true, size_t BUCKET_OFFSET=8>
class khpp {
    // BUCKET_OFFSET is
    using index_type = u64;
    index_type n_buckets_, size_, n_used_;

    uint32_t *flags_;
    using pair_type = typename khel_t<khkey_t, khval_t, is_map>::Type;
    pair_type *pairs_;
    hash_equal he;
    hash_func  hf;
    std::vector<tthread::fast_mutex> locks;
    tthread::fast_mutex global_lock_;

public:

    static const size_t LOCK_BUCKET_SIZE = 1ull << BUCKET_OFFSET;
    static constexpr double HASH_UPPER = 0.77;

    size_t flag_size(size_t m) {
        return m < 16 ? 1 : m >> 4;
    }
    size_t nlocks(size_t num_buckets) const {
        return num_buckets < (1 << BUCKET_OFFSET) ? 1: num_buckets >> BUCKET_OFFSET;
    }
    auto upper_bound() const {
        return (index_type)(n_buckets_ * HASH_UPPER + 0.5);
    }
    auto is_empty(index_type i) const {
        return (flags_[i>>4]>>((i&0xfU)<<1))&2;
    }
    auto is_del(index_type i) const {
        return (flags_[i>>4]>>((i&0xfU)<<1))&1;
    }
    auto is_either(index_type i) const {
        return (flags_[i>>4]>>((i&0xfU)<<1))&3;
    }
    auto exists(index_type i) const {
        return (((flags_[i>>4]>>((i&0xfU)<<1)))&3) == 0;
    }
    u64 estimate_memory() {
        return n_buckets_ * sizeof(pair_type) + flag_size(n_buckets_) * sizeof(uint32_t);
    }

    using key_type = khkey_t;
    using val_type = khval_t;
    khpp() {
        std::memset(this, 0, sizeof(*this));
    }
    khpp(size_t size): khpp() {
        resize(size);
    }
    void clear() {
        global_lock_.lock();
        if constexpr(!pair_type::trivially_destructible()) {
            for(index_type i(0); i < n_buckets_; ++i) {
                if(exists(i)) {
                    pairs_[i].~pair_type();
                }
            }
        }
        std::memset(flags_, 0xaa, flag_size(n_buckets_) * sizeof(uint32_t));
        size_ = n_used_ = 0;
        global_lock_.unlock();
    }
    ~khpp() {
        global_lock_.lock();
        free(pairs_);
        free(flags_);
        global_lock_.unlock();
    }
	int resize(index_type new_n_buckets)
	{ /* This function uses 0.25*n_buckets bytes of working space instead of [sizeof(key_t+val_t)+.25]*n_buckets. */
        global_lock_.lock();
		uint32_t *new_flags = 0;
		u64 j = 1;
		{
			kroundup64(new_n_buckets);
			if (new_n_buckets < 4) new_n_buckets = 4;
			if (size_ >= (u64)(new_n_buckets * HASH_UPPER + 0.5)) j = 0;	/* requested size is too small */
			else { /* hash table size to be changed (shrink or expand); rehash */
				new_flags = (uint32_t*)malloc(flag_size(new_n_buckets) * sizeof(uint32_t));
				if (!new_flags) {
                    global_lock_.unlock();
                    return -1;
                }
				memset(new_flags, 0xaa, flag_size(new_n_buckets) * sizeof(uint32_t));
				if (n_buckets_ < new_n_buckets) {	/* expand */
                    pair_type *new_pairs(static_cast<pair_type *>(realloc((void *)pairs_, new_n_buckets * sizeof(pair_type))));
					if (!new_pairs) { free(new_flags); global_lock_.unlock(); return -1; }
				} /* otherwise shrink */
			}
		}
		if (j) { /* rehashing is needed */
			for (j = 0; j < n_buckets_; ++j) {
                if(exists(j)) {
                    pair_type &pair(pairs_[j]);
					const u64 new_mask(new_n_buckets - 1);
					__ac_set_isdel_true(flags_, j);
					while (1) { /* kick-out process; sort of like in Cuckoo hashing */
						u64 k, i, step = 0;
						k = hf(pair.first);
						i = k & new_mask;
						while (!__ac_isempty(new_flags, i)) i = (i + (++step)) & new_mask;
						__ac_set_isempty_false(new_flags, i);
						if (i < n_buckets_ && __ac_iseither(flags_, i) == 0) { /* kick out the existing element */
							{
                                auto tmp(std::move(pair.first));
                                pair.first = std::move(pairs_[i].first);
                                pairs_[i].first = std::move(tmp);
                            }
							if constexpr(is_map) {
                                auto tmp(std::move(pair.second));
                                pair.second = std::move(pairs_[i].second);
                                pairs_[i].second = std::move(tmp);
                            }
							__ac_set_isdel_true(flags_, i); /* mark it as deleted in the old hash table */
						} else { /* write the element and jump out of the loop */
                            pairs_[i] = std::move(pair);
							break;
						}
					}
				}
			}
			if (n_buckets_ > new_n_buckets) { /* shrink the hash table */
                pairs_ = (pair_type *)realloc((void *)pairs_, new_n_buckets * sizeof(pair_type));
			}
			free(flags_); /* free the working space */
			flags_ = new_flags;
			n_buckets_ = new_n_buckets;
			n_used_ = size_;
		}
        global_lock_.unlock();
        locks.resize(nlocks(n_buckets_));
		return 0;
	}
    template<typename... Args>
    pair_type &emplace(khkey_t &key, int *ret)
    {
        index_type x;
        if (n_occupied >= upper_bound) { /* update the hash table */
            if (n_buckets > (size<<1)) {
                if (resize(h, n_buckets - 1) < 0) { /* clear "deleted" elements */
                    *ret = -1; return n_buckets;
                }
            } else if (resize(h, n_buckets + 1) < 0) { /* expand the hash table */
                *ret = -1; return n_buckets;
            }
        } /* TODO: to implement automatically shrinking; resize() already support shrinking */
        {
            index_type k, i, site, last, mask = n_buckets - 1, step = 0;
            x = site = n_buckets; k = __hash_func(key); i = k & mask;
            if (__ac_isempty(flags, i)) x = i; /* for speed up */
            else {
                last = i;
                while (!__ac_isempty(flags, i) && (__ac_isdel(flags, i) || !__hash_equal(keys[i], key))) {
                    if (__ac_isdel(flags, i)) site = i;
                    i = (i + (++step)) & mask;
                    if (i == last) { x = site; break; }
                }
                if (x == n_buckets) {
                    if (__ac_isempty(flags, i) && site != n_buckets) x = site;
                    else x = i;
                }
            }
        }
        std::unique_lock<std::shared_mutex>(m);
        if (__ac_isempty(flags, x)) { /* not present at all */
            keys[x] = key;
            __ac_set_isboth_false(flags, x);
            ++size; ++n_occupied;
            *ret = 1;
        } else if (__ac_isdel(flags, x)) { /* deleted */
            keys[x] = key;
            __ac_set_isboth_false(flags, x);
            ++size;
            *ret = 2;
        } else *ret = 0; /* Don't touch keys[x] if present and not deleted */
        return x;
    }
};

} // namespace kh

#endif
