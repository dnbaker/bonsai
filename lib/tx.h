#ifndef TX_H__
#define TX_H__

#include "lib/feature_min.h"
#include "lib/util.h"

namespace emp {


class Taxonomy {
    khash_t(p) *tax_map_;
    khash_t(name) *name_map_;
    uint64_t synthetic_count_, synthetic_ceiling_;
public:
    // Textual constructor
    Taxonomy(const char *taxnodes_path, const char *name_path):
        tax_map_(khash_load<khash_t(p)>(taxnodes_path)),
        name_map_(build_name_hash(name_path)),
        synthetic_count_(0),
        synthetic_ceiling(-1)
    {
    }
    // constructor
    ~Taxonomy() {
        kh_destroy(p, tax_map_);
        destroy_name_hash(name_map_);
    }
    void write(const char *fn);
    void add_node_impl(const char *node_name, const unsigned node_id, const unsigned parent);
    void add_node(const char *node_name, const unsigned parent);
};

template<typename T>
class EMSetTree {
    std::vector<T *> core_;
};

} // namespace emp

#endif // #ifndef TX_H__
