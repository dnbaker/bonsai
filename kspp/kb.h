#ifndef KBTREE_WRAPPER_H__
#define KBTREE_WRAPPER_H__
#include "klib/kbtree.h"
#include <cstring>
#include <stdexcept>
#include <string>
#include <cassert>

#ifndef unlikely
#define unlikely(x) __builtin_expect((x), 0)
#endif
#ifndef likely
#define likely(x) __builtin_expect((x), 1)
#endif


namespace kb {

struct DefaultCmp {
    template<typename T>
    int operator()(const T &a, const T &b) const {
        return (a < b) ? -1: b < a ? 1: 0;
    }
};

template<typename KeyType, typename Cmp=DefaultCmp, size_t MAX_DEPTH_PARAM=64>
class KBTree {
public:
    using key_t = KeyType;
    using pointer = key_t *;
    using const_pointer = const key_t *;
    static constexpr size_t MAX_DEPTH = MAX_DEPTH_PARAM;

    struct node_t {int32_t is_internal:1, n:31;};
    struct pos_t  {node_t *x; int i;};
    struct iter_t {
        pos_t stack[MAX_DEPTH], *p;

        iter_t(KBTree &ref): p(nullptr) {
            if(unlikely(ref.n_keys == 0)) throw std::runtime_error("Could not initialize iterator over empty sequence.");
            p = stack;
            p->x = ref.root; p->i = 0;
            while(p->x->is_internal && ref.ptr(p->x)[0] != 0) {
                node_t *x = p->x;
                (++p)->x = ref.ptr(x)[0];
                p->i = 0;
            }
        }
        key_t &key() {
            return reinterpret_cast<key_t *>(reinterpret_cast<char *>(p->x) + 4)[p->i];
        }
        const key_t &key() const {
            return reinterpret_cast<key_t *>(reinterpret_cast<char *>(p->x) + 4)[p->i];
        }
        const key_t &const_key() const {
#if !NDEBUG
            auto ptr = reinterpret_cast<key_t *>(reinterpret_cast<char *>(p->x) + 4) + p->i;
            if(ptr == nullptr) {
                std::fprintf(stderr, "Warning: Null key at position\n");
            }
#endif
            return reinterpret_cast<key_t *>(reinterpret_cast<char *>(p->x) + 4)[p->i];
        }
        bool valid() const {return p >= stack;}
    };

    int t, n;
    int off_key, off_ptr, ilen, elen;
    node_t *root;
    int n_keys, n_nodes;
    KBTree(size_t size=KB_DEFAULT_SIZE):
        t(((size - 4 - sizeof(void *)) / (sizeof(void *) + sizeof(key_t)) + 1) >>1), n((t<<1) - 1),
        off_ptr(4 + n * sizeof(key_t)), ilen((4 + sizeof(void *) + n * (sizeof(void *) + sizeof(key_t)) + 3) >> 2 << 2),
        root(t >= 2 ? static_cast<node_t *>(std::calloc(1, ilen)): nullptr),
        off_key(0),
        elen((off_ptr + 3) >> 2 << 2), n_keys(0), n_nodes(1)
    {
        if(!root) throw std::invalid_argument(std::string("t must be >= 2. t: ") + std::to_string(t));
    };
    node_t **ptr(node_t *x) {
        return reinterpret_cast<node_t **>(reinterpret_cast<char *>(this) + off_ptr);
    }
    const node_t **ptr(node_t *x) const {
        return reinterpret_cast<const node_t **>(reinterpret_cast<const char *>(this) + off_ptr);
    }
    int proot(std::FILE *fp=stderr) const {return std::fprintf(fp, "root: %p\n", root);}
    ~KBTree() {
        std::fprintf(stderr, "Calling destructor\n");
#if 0
        top = stack = (kbnode_t**)calloc(max, sizeof(kbnode_t*));
            *top++ = (b)->root;
            while (top != stack) {
                x = *--top;
                if (x->is_internal == 0) { free(x); continue; }
                for (i = 0; i <= x->n; ++i)
                    if (__KB_PTR(b, x)[i]) {
                        if (top - stack == max) {
                            max <<= 1;
                            stack = (kbnode_t**)realloc(stack, max * sizeof(kbnode_t*));
                            top = stack + (max>>1);
                        }
                        *top++ = __KB_PTR(b, x)[i];
                    }
                free(x);
            }
        }
        free(stack);
#endif
        int i, max = 8;
        node_t *x, **top, **stack = 0;
        std::fprintf(stderr, "Allocating workspace\n");
        proot(stderr);
        top = stack = static_cast<node_t **>(std::calloc(max, sizeof(node_t *)));
        assert(stack); assert(top);
        std::fprintf(stderr, "Further down. top ptr: %p, root %p\n", static_cast<void *>(top), static_cast<void *>(root));
        *top++ = root;
        std::fprintf(stderr, "root ptr: %p\n", static_cast<void *>(root));
        std::fflush(stderr);
        while(top != stack) {
            x = *--top;
            assert(x);
            if(x->is_internal == 0) {std::free(x); continue;}
            for(i = 0; i <= n; ++i) {
                if (ptr(x)[i]) {
                    if(top - stack == max) {
                        max <<= 1;
                        stack = static_cast<node_t **>(std::realloc(stack, max * sizeof(node_t *)));
                    }
                    *top++ = ptr(x)[i];
                }
                std::free(x);
            }
            std::free(stack);
        }
    }
    int cmp(const KeyType &a, const KeyType &b) const {return Cmp()(a, b);}
    static key_t *key(node_t *node) {
        return reinterpret_cast<key_t *>(reinterpret_cast<char *>(node) + 4);
    }
    static const key_t *key(const node_t *node) {
        return reinterpret_cast<const key_t *>(reinterpret_cast<const char *>(node) + 4);
    }
    static int get_aux(const node_t * __restrict x, const key_t * __restrict k, int *r) {
        int tr, *rr, begin = 0, end = x->n;
        if (x->n == 0) return -1;
        rr = r ? r: &tr;
        while(begin < end) {
            int mid = (begin + end) >> 1;
            if(Cmp()(key(x)[mid], *k) < 0) begin = mid + 1;
            else end = mid;
        }
        if(begin == x->n) { *rr = 1; return x->n - 1;}
        begin -= (*rr = Cmp()(*k, key(x)[begin])) < 0;
        return begin;
    }
    key_t *get(const key_t * __restrict k) {
        int i, r = 0;
        node_t *x = root;
        while (x) {
            i = get_aux(x, k, &r);
            if(i >= 0 && r == 0) return &key(x)[i];
            if (x->is_internal == 0) return nullptr;
            x = ptr(x)[i + 1];
        }
        return nullptr;
    }
    const key_t *get(const key_t * __restrict k) const {
        int i, r = 0;
        node_t *x = root;
        while (x) {
            i = get_aux(x, k, &r);
            if(i >= 0 && r == 0) return &key(x)[i];
            if (x->is_internal == 0) return nullptr;
            x = ptr(x)[i + 1];
        }
        return nullptr;
    }
    key_t *get(key_t k) {return get(&k);}
    const key_t *get(key_t k) const {return get(&k);}
    auto size() const {return n_keys;}
    void interval(const key_t * __restrict k, key_t **lower, key_t **upper) {
        int i, r = 0;
        node_t *x = root;
        *lower = *upper = 0;
        while (x) {
            i = get_aux(x, k, &r);
            if (i >= 0 && r == 0) {
                *lower = *upper = &key(x)[i];
                return;
            }
            if (i >= 0) *lower = &__KB_KEY(key_t, x)[i];
            if (i < x->n - 1) *upper = &key(x)[i + 1];
            if (x->is_internal == 0) return;
            x = ptr(x)[i + 1];
        }
    }
    void interval(const key_t k, key_t **lower, key_t **upper) {
        interval(&k, lower, upper);
    }
    void split(node_t *x, int i, node_t *y) {
        node_t *z = static_cast<node_t *>(std::calloc(1, y->is_internal ? ilen: elen));
        ++n_nodes;
        z->is_internal = y->is_internal;
        z->n = this->t - 1;
        std::memcpy(key(z), key(y) + this->t, sizeof(key_t) * (this->t - 1));
        if(y->is_internal) std::memcpy(ptr(z), ptr(y) + this->t, sizeof(void *) * this->t);
        y->n = this->t - 1;
        std::memmove(ptr(x) + i + 2, ptr(x) + i + 1, sizeof(void *) * (x->n - i));
        ptr(x)[i + 1] = z;
        std::memmove(key(x) + i + 1, key(x) + i, sizeof(key_t) * (x->n - i));
        key(x)[i] = key(y)[this->t - 1];
        ++x->n;
    }
    key_t *put_aux(node_t *x, const key_t * __restrict k) {
        int i = x->n - 1;
        key_t *ret;
        if (x->is_internal == 0) {
            if((i = get_aux(x, k, 0)) != x->n - 1)
                std::memmove(key(x) + i + 2, key(x) + i + 1, (x->n - i - 1) * sizeof(key_t));
            ret = &key(x)[i + 1];
            *ret = *k;
            ++x->n;
        } else {
            i = get_aux(x, k, 0) + 1;
            if (ptr(x)[i]->n == 2 * this->t - 1) {
                split(x, i, ptr(x)[i]);
                i += (Cmp()(*k, key(x)[i]) > 0);
            }
            ret = put_aux(ptr(x)[i], k);
        }
        return ret;
    }
    key_t *put(const key_t k) {
        return put(&k);
    }
    key_t *put(const key_t * __restrict k) {
        node_t *r, *s;
        ++this->n_keys;
        r = this->root;
        if (r->n == 2 * this->t - 1) {
            ++this->n_nodes;
            this->root = s = static_cast<node_t*>(calloc(1, this->ilen));
            s->is_internal = 1; s->n = 0;
            ptr(s)[0] = r;
            split(s, 0, r);
            r = s;
        }
        return put_aux(r, k);
    }
    int itr_get(const key_t * __restrict k, iter_t &itr) {
        int i, r = 0;
        itr.p = itr.stack;
        itr.p->x = root; itr.p->i = 0;
        while(itr.p->x) {
            i = get_aux(itr.p->x, k, &r);
            if (i >= 0 && r == 0) return 0;
            if (itr.p->x->is_internal == 0) return -1;
            itr.p[1].x = ptr(itr.p->x)[i + 1];
            itr.p[1].i = i;
            ++itr.p;
        }
    }
    int itr_next(iter_t &itr) {
        if(itr.p < itr.stack) return 0;
        for(;;) {
            ++itr.p->i;
            while (itr.p->x && itr.p->i <= itr.p->x->n) {
                itr.p[1].i = 0;
                itr.p[1].x = itr.p->x->is_internal? ptr(itr.p->x)[itr.p->i] : 0;
                ++itr.p;
            }
            --itr.p;
            if (itr.p < itr.stack) return 0;
            if (itr.p->x && itr.p->i < itr.p->x->n) return 1;
        }
    }
    int itr_next(iter_t *itr) {return itr_next(*itr);}
    template<typename Func>
    void for_each(const Func &func) {
        proot(stderr);
        for(iter_t it(*this); it.valid(); func(it.key()), itr_next(it));
    }
    template<typename Func>
    void for_each(const Func &func) const {
        proot(stderr);
        for(iter_t it(*this); it.valid(); func(it.const_key()), itr_next(it));
    }
    // TODO: Delete
};


} // namespace kb

#endif
