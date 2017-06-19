#ifndef _UNION_FIND__
#define _UNION_FIND__
#include <vector>
#include <cstdio>

#if (defined(__GNUC__) && __GNUC__) || (defined(__clang__) && __clang__)
#define PACKED __attribute__((packed))
#else
#define PACKED
#endif

namespace uf { namespace fds {

template<typename T, typename size_type=std::uint32_t>
class DisjointSetVector {
    struct ufnode_t {
        T          e_; // element
        size_type  p_; // parent
        unsigned   r_;
        template<typename... Args>
        ufnode_t(size_type index, Args&&... args):
            e_(std::forward<T>(args)...), p_(index), r_{0u} {
        }
    } PACKED;
#if !NDEBUG
public:
#endif
    std::vector<ufnode_t>  v_;
    size_type           n_[6];
#if NDEBUG
public:
#endif
    size_type find(ufnode_t &node) {
        const size_type index(&node - v_.data());
        return find(&node - v_.data());
    }
    size_type find(size_type index) {
        if(v_[index].p_ == index) v_[index].p_ = find(v_[v_[index].p_]);
        return v_[index].p_;
    }

    template<typename... Args>
    size_type emplace(Args&&... args) {
        v_.emplace_back(args...);
        ++n_[0];
        return v_.size() - 1;
    }

    auto size()   const {return v_.size();}
    auto begin()        {return v_.begin();}
    auto end()          {return v_.end();}
    auto cbegin() const {return v_.cbegin();}
    auto cend()   const {return v_.cend();}
    DisjointSetVector():  n_{0} {}

    void perform_union(size_type ai, size_type bi) {
        if(ai == bi) return;
        if(v_[ai].r_ < v_[bi].r_) v_[ai].p_ = bi;
        if(v_[ai].r_ > v_[bi].r_) v_[bi].p_ = ai;
        else v_[ai].p_ = bi, --n_[v_[ai].r_], ++n_[++v_[bi].r_];
    }

    void perform_union(ufnode_t &a, ufnode_t &b) {
        return perform_union(find(a), find(b));
    }

    bool same_set(ufnode_t &a, ufnode_t &b) {return find(a) == find(b);}
};

}} // namespace emp::fds

#endif // _UNION_FIND__
