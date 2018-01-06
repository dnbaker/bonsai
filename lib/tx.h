#ifndef _TX_H__
#define _TX_H__
#include <stdexcept>
#include <thread>
#include <unordered_map>
#include <functional>
#include "klib/kthread.h"
#include "lib/feature_min.h"
#include "lib/counter.h"
#include "lib/bits.h"


template<typename T>
size_t get_n_occ(T *hash) {
    size_t ret(0);
    for(khiter_t ki(0); ki != kh_end(hash); ++ki) ret += !!kh_exist(hash, ki);
    return ret;
}


namespace emp {


using namespace std::literals;

/*
1. Get taxids you need.
 */
class TaxonomyReformation {
protected:
    khash_t(p)           *pmap_;
    khash_t(p)           *rpmap_;
    khash_t(name)    *name_map_;
    std::vector<tax_t> old_ids_;
    std::unordered_map<std::string, tax_t> path_map;
    uint64_t        counter_:63;
    uint64_t          filled_:1;
public:
    template<typename StrCon>
    TaxonomyReformation(const char *name_path, const StrCon &paths):
        pmap_{kh_init(p)}, rpmap_{kh_init(p)},
        name_map_{build_name_hash(name_path)},
        old_ids_{0, 1},
        counter_(1), filled_(false)
    {
        // Add the root of the tree.
        int khr;
        khint_t ki{kh_put(p, pmap_, tax_t(counter_), &khr)};
        kh_val(pmap_, ki) = 0;
        ki = kh_put(p, rpmap_, 0, &khr);
        kh_val(rpmap_, ki) = 1;
        fill_path_map(paths);
        count::Counter<tax_t> ctr;
        for(const auto &pair: path_map) ctr.add(pair.second);
    }

    template<typename T>
    void fill_path_map(const T &container) {
        for(const std::string &path: container) {
            const tax_t id(get_taxid(get_str(path), name_map_));
            if(id == tax_t(-1)) {
                LOG_WARNING("Tax id not found in path %s. Skipping. This can be fixed by augmenting the name dictionary file.\n", path);
                continue;
            }
            path_map[path] = id;
        }
    }
    
    void add_genome(const char *path) {
        assert(name_map_);
        const tax_t id(get_taxid(path, name_map_));
        if(id == tax_t(-1)) {
            LOG_WARNING("Tax id not found in path %s. Skipping. This can be fixed by augmenting the name dictionary file.\n", path);
            return;
        }
        khint_t ki;
        if((ki = kh_get(p, pmap_, id)) == kh_end(pmap_)) {
            // Add the node, insert its parent (I will also need access to the NCBI taxonomy for that.)
        } else {
            // This node is now 
        }
        if(parent(id) == tax_t(-1)) LOG_EXIT("Error in taxonomy hash map.");
        assert(!"NotImplemented.");
    }
    template<typename StrType> void add_genome(const StrType &path) {add_genome(path.data());}

    // Convenience functions
    tax_t parent(tax_t child) const {
        const auto ind(kh_get(p, pmap_, child));
        return ind == kh_end(pmap_) ? tax_t(-1): kh_val(pmap_, ind);
    }
    ~TaxonomyReformation() {
        khash_destroy(name_map_);
        khash_destroy(pmap_);
        khash_destroy(rpmap_);
    }
};
using concise_tax_t = TaxonomyReformation;



} // namespace emp

#endif // #ifndef _TX_H__
