#ifndef TX_H__
#define TX_H__
#include <stdexcept>

#include "lib/feature_min.h"
#include "lib/util.h"

namespace emp {

class Taxonomy {
    khash_t(p) *tax_map_;
    khash_t(name) *name_map_;
    uint64_t n_syn_;
    uint64_t ceil_;
public:
    // Textual constructor
    Taxonomy(const char *taxnodes_path, const char *name_path, unsigned ceil=0):
        tax_map_(khash_load<khash_t(p)>(taxnodes_path)),
        name_map_(build_name_hash(name_path)),
        n_syn_(0),
        ceil_(ceil ? ceil: tax_map_->n_buckets << 1)
    {
    }
    // Binary constructor
    Taxonomy(const char *path, unsigned ceil=0);
    ~Taxonomy() {
        kh_destroy(p, tax_map_);
        destroy_name_hash(name_map_);
    }
    void write(const char *fn) const;
    void add_node_impl(const char *node_name, const unsigned node_id, const unsigned parent);
    void add_node(const char *node_name, const unsigned parent);
    uint64_t get_syn_count() const {return n_syn_;}
    uint64_t get_syn_ceil() const {return ceil_;}
    int can_add() { return n_syn_ + 1 <= ceil_;}
};

template<typename T>
class EMSetTree {
    std::vector<T *> core_;
};

} // namespace emp

#endif // #ifndef TX_H__
