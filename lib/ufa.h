#ifndef _UNION_FIND_ADAPTER__
#define _UNION_FIND_ADAPTER__

#if (defined(__GNUC__) && __GNUC__) || (defined(__clang__) && __clang__)
#define PACKED __attribute__((packed))
#else
#define PACKED
#endif

namespace emp { namespace ufa {

// For composition
template<typename T, typename size_type=std::uint8_t>
struct uf_adapter: public T {
    size_type      r_;
    uf_adapter    *p_;
    template<typename... Args>
    uf_adapter(Args&&... args):
        T(std::forward<Args...>(args)...), r_{0}, p_{this} {}
} PACKED;

template<typename T>
T *find(T *node) {
    if(node->p_ != node) node->p_ = find(node->p_);
    return node->p_;
}

template<typename T>
T *find(T &node) {return find(std::addressof(node));}

template<typename T>
void perform_union(T *a, T *b) {
    if((a = find(a)) == (b = find(b))) return;
    if     (a->r_ < b->r_) b->p_ = a;
    else if(b->r_ < b->r_) a->p_ = b;
    else          b->p_ = a, ++a->r_;
}

template<typename T>
void perform_union(T &a, T &b) {perform_union(std::addressof(a), std::addressof(b));}

}} // namespace emp::ufa

#endif // _UNION_FIND_ADAPTER__
