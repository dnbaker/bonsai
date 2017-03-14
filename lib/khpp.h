#ifndef KHPP_H__
#define KHPP_H__
#include "lib/khash64.h"
#include "lib/hash.h"

namespace emp {
namespace kh {

template<typename khkey_t, typename khval_t, typename hash_func=std::hash<khkey_t>, class hash_equal=std::equal_to<khkey_t>, bool is_map=true>
struct khpp_t {
    using index_type = std::uint64_t;
    index_type n_buckets, size, n_occupied, upper_bound;
    khint32_t *flags;
    khkey_t *keys;
    khval_t *vals;
    mutable std::shared_mutex m;
    hash_equal he;
    hash_func hf;
    static constexpr double HASH_UPPER = 0.77;
    khpp_t(): n_buckets(0), size(0), n_occupied(0), upper_bound(0), flags(0), keys(0), vals(0) {
    }
    ~khpp_t() {
        free(flags);
        free(vals);
        free(keys);
    }
    void clear() {
        std::unique_lock<std::shared_mutex>(m);
        if (flags) {
            memset(flags, 0xaa, __ac_fsize(n_buckets) * sizeof(khint32_t));
            size = n_occupied = 0;
        }
    }
    index_type iget(khkey_t &key)
    {
        std::shared_lock<std::shared_mutex>(m);
        if (n_buckets) {
            index_type k, i, last, mask, step = 0;
            mask = n_buckets - 1;
            k = hf(key); i = k & mask;
            last = i;
            while (!__ac_isempty(flags, i) && (__ac_isdel(flags, i) || !he(keys[i], key))) {
                i = (i + (++step)) & mask;
                if (i == last) return n_buckets;
            }
            return __ac_iseither(flags, i)? n_buckets : i;
        }
        return 0;
    }
    index_type iput(khkey_t &key, int *ret)
    {
        index_type x;
        if (n_occupied >= upper_bound) { /* update the hash table */
            if (n_buckets > (size<<1)) {
                if (resize(n_buckets - 1) < 0) { /* clear "deleted" elements */
                    *ret = -1; return n_buckets;
                }
            } else if (resize(n_buckets + 1) < 0) { /* expand the hash table */
                *ret = -1; return n_buckets;
            }
        } /* TODO: to implement automatically shrinking; resize() already support shrinking */
        std::unique_lock<std::shared_mutex>(m);
        {
            index_type k, i, site, last, mask = n_buckets - 1, step = 0;
            x = site = n_buckets; k = hf(key); i = k & mask;
            if (__ac_isempty(flags, i)) x = i; /* for speed up */
            else {
                last = i;
                while (!__ac_isempty(flags, i) && (__ac_isdel(flags, i) || !he(keys[i], key))) {
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
    index_type nb() const {return n_buckets;}
    khval_t &operator[](khkey_t &key) {
        std::shared_lock<std::shared_mutex>(m);
        index_type ki(iget(key));
        int khr;
        if(ki == nb()) ki = iput(key, &khr);
        return vals[ki];
    }
    void del(index_type x)
    {
        std::unique_lock<std::shared_mutex>(m);
        if (x != n_buckets && !__ac_iseither(flags, x)) {
            __ac_set_isdel_true(flags, x);
            --size;
        }
    }
    int resize(index_type new_n_buckets)
    { /* This function uses 0.25*n_buckets bytes of working space instead of [sizeof(key_t+val_t)+.25]*n_buckets. */
        khint32_t *new_flags = 0;
        index_type j = 1;
        {
            kroundup64(new_n_buckets);
            if (new_n_buckets < 4) new_n_buckets = 4;
            if (size >= (index_type)(new_n_buckets * HASH_UPPER + 0.5)) j = 0;    /* requested size is too small */
            else { /* hash table size to be changed (shrink or expand); rehash */
                std::unique_lock<std::shared_mutex>(m);
                new_flags = (khint32_t*)kmalloc(__ac_fsize(new_n_buckets) * sizeof(khint32_t));
                if (!new_flags) return -1;
                memset(new_flags, 0xaa, __ac_fsize(new_n_buckets) * sizeof(khint32_t));
                if (n_buckets < new_n_buckets) {    /* expand */
                    khkey_t *new_keys = (khkey_t*)krealloc((void *)keys, new_n_buckets * sizeof(khkey_t));
                    if (!new_keys) { kfree(new_flags); return -1; }
                    keys = new_keys;
                    if (is_map) {
                        khval_t *new_vals = (khval_t*)krealloc((void *)vals, new_n_buckets * sizeof(khval_t));
                        if (!new_vals) { kfree(new_flags); return -1; }
                        vals = new_vals;
                    }
                } /* otherwise shrink */
            }
        }
        if (j) { /* rehashing is needed */
            std::unique_lock<std::shared_mutex>(m);
            for (j = 0; j != n_buckets; ++j) {
                if (__ac_iseither(flags, j) == 0) {
                    khkey_t key = keys[j];
                    khval_t val;
                    index_type new_mask;
                    new_mask = new_n_buckets - 1;
                    if (is_map) val = vals[j];
                    __ac_set_isdel_true(flags, j);
                    while (1) { /* kick-out process; sort of like in Cuckoo hashing */
                        index_type k, i, step = 0;
                        k = hf(key);
                        i = k & new_mask;
                        while (!__ac_isempty(new_flags, i)) i = (i + (++step)) & new_mask;
                        __ac_set_isempty_false(new_flags, i);
                        if (i < n_buckets && __ac_iseither(flags, i) == 0) { /* kick out the existing element */
                            { khkey_t tmp = keys[i]; keys[i] = key; key = tmp; }
                            if (is_map) { khval_t tmp = vals[i]; vals[i] = val; val = tmp; }
                            __ac_set_isdel_true(flags, i); /* mark it as deleted in the old hash table */
                        } else { /* write the element and jump out of the loop */
                            keys[i] = key;
                            if (is_map) vals[i] = val;
                            break;
                        }
                    }
                }
            }
            if (n_buckets > new_n_buckets) { /* shrink the hash table */
                keys = (khkey_t*)krealloc((void *)keys, new_n_buckets * sizeof(khkey_t));
                if (is_map) vals = (khval_t*)krealloc((void *)vals, new_n_buckets * sizeof(khval_t));
            }
            kfree(flags); /* free the working space */
            flags = new_flags;
            n_buckets = new_n_buckets;
            n_occupied = size;
            upper_bound = (index_type)(n_buckets * HASH_UPPER + 0.5);
        }
        return 0;
    }
    bool try_set(const index_type ki, const khval_t &val)
    {
        std::shared_lock<std::shared_mutex>(m);
        return __sync_bool_compare_and_swap(vals + ki, vals[ki], val);
    }
    void set(khkey_t &key, const khval_t &val)
    {
        khiter_t ki;
        int khr;
        if((ki = iget(key)) == nb()) ki = iput(key, &khr);
        while(!try_set(ki, val));
    }
    template<typename T>
    void func_set(T func, khkey_t &key, const khval_t &val) {
        std::shared_lock<std::shared_mutex> lock(m);
        func_set_impl(func, key, val);
    }
    template<typename T>
    void func_set_impl(T func, khkey_t &key, const khval_t &val) {
        std::shared_lock<std::shared_mutex> lock(m);
        func_set_impl(func, key, val);
    }

};

}
}
#endif
