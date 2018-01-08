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
#include "lib/linearset.h"


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
2. Create the map over and the map back.
3. Move the new id/file path mapping into workflow for logging.
 */
class TaxonomyReformation {
protected:
    khash_t(p)           *pmap_;
    khash_t(name)    *name_map_;
    std::vector<tax_t> old_ids_;
    std::unordered_map<tax_t, std::vector<std::string>> path_map;
    std::unordered_map<tax_t, std::vector<std::string>> newid_path_map;
    uint64_t        counter_:62;
    uint64_t          filled_:1;
    uint64_t  panic_on_undef_:1;
public:
    template<typename StrCon>
    TaxonomyReformation(const char *name_path, const StrCon &paths,
                        const khash_t(p) *old_tax, bool panic_on_undef=false):
        pmap_{kh_init(p)},
        name_map_{build_name_hash(name_path)},
        old_ids_{0, 1},
        counter_(1), filled_(false), panic_on_undef_(panic_on_undef)
    {
        khash_t(p) *ct(kh_init(p));
        kh_resize(p, ct, kh_size(old_tax));
        std::memcpy(ct->keys, old_tax->keys, sizeof(tax_t) * old_tax->n_buckets);
        std::memcpy(ct->vals, old_tax->vals, sizeof(tax_t) * old_tax->n_buckets);
        std::memcpy(ct->flags, old_tax->flags, __ac_fsize(old_tax->n_buckets));
#define CP(x) ct->x = old_tax->x
        CP(n_buckets); CP(size); CP(n_occupied); CP(upper_bound);
#undef CP
        // Add the root of the tree.
        int khr;
        khint_t ki(kh_put(p, pmap_, tax_t(counter_), &khr));
        kh_val(pmap_, ki) = 0;
        fill_path_map(paths);
#if CHUNKY_BUT_FUNKY
        std::unordered_set<tax_t> funkynodes; // Nodes with more than one genome which we'll have to break up.
        for(const auto &pair: path_map)
            if(pair.second.size() > 1)
                funkynodes.insert(pair.first);
#endif
        // Arbitrary seeding scheme which still will vary run to run
        std::mt19937 mt(std::hash<size_t>()(kh_size(old_tax) * path_map.size()));
        std::unordered_map<tax_t, tax_t> old_to_new;
        for(auto &pair: path_map) {
            if(pair.second.size() > 1) {
                for(auto &path: pair.second) {
                    tax_t tmp(mt());
                    while(kh_get(p, ct, tmp) != kh_end(ct)) tmp = mt();
                    path_map[tmp] = {path};
                    ki = kh_put(p, ct, tmp, &khr);
                    kh_val(ct, ki) = pair.first;
                    newid_path_map[tmp] = {path};
                }
                std::vector<std::string> vec;
                std::swap(vec, pair.second);
            }
        }
        std::vector<tax_t> insertion_order;
        insertion_order.reserve(path_map.size());
        for(const auto &[tax, _]: path_map) insertion_order.push_back(tax);
        assert(insertion_order.size() == std::unordered_set<tax_t>(insertion_order.begin(), insertion_order.end()).size() &&
               "Some tax ids are being repeated in this array.");
        pdqsort(std::begin(insertion_order), std::end(insertion_order),
                [old_tax] (const tax_t a, const tax_t b) {
            return node_depth(old_tax, a) < node_depth(old_tax, b);
        });
        for(const auto tax: insertion_order) {
            old_to_new[tax] = old_ids_.size();
            old_ids_.emplace_back(tax);
            ki = kh_put(p, pmap_, static_cast<tax_t>(old_ids_.size()), &khr);
            kh_val(pmap_, ki) = old_to_new.at(kh_val(ct, kh_get(p, ct, tax)));
            // TODO: remove \.at check when we're certain it's working.
        }
        khash_destroy(ct);
        // Now convert name_map from old to new.
        for(khint_t ki(0); ki < kh_size(name_map_); ++ki)
            if(kh_exist(name_map_, ki))
                kh_val(name_map_, ki) = \
                    old_to_new.at(kh_val(name_map_, ki));
                    // TODO: remove \.at check when we're certain it's working.

        LOG_DEBUG("Paths to genomes with new subtax elements:\n\n\n%s", newtaxprintf().data());
        
    }

    ks::string newtaxprintf() {
        ks::string ret("#New ID\tGenome path (NEW FIRST)\n");
        for(const auto &pair: newid_path_map) {
            ret.sprintf("%u\t%s\n", pair.first, pair.second.data());
        }
        ret += "#(Original taxonomy)\n";
        for(const auto &pair: path_map) {
            ret.sprintf("%u\t%s\n", pair.first, pair.second.data());
        }
        return ret;
    }
    void fnewtaxprintf(std::FILE *fp) {
        auto taxstr(newtaxprintf());
        std::fwrite(taxstr.data(), 1, taxstr.size(), fp);
    }

    template<typename T>
    void fill_path_map(const T &container) {
        for(const auto &path: container) {
            const tax_t id(get_taxid(get_str(path), name_map_));
            if(id == tax_t(-1)) {
                if(panic_on_undef_)
                    throw std::runtime_error(ks::sprintf("Tax id not found in path %s. Skipping. This can be fixed by augmenting the name dictionary file.\n", get_str(path)).data());
                else
                    LOG_WARNING("Tax id not found in path %s. Skipping. This can be fixed by augmenting the name dictionary file.\n", get_str(path));
                continue;
            }
            path_map[id].push_back(path);
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
    template<typename StrType> void add_genome(const StrType &path) {add_genome(get_str(path));}

    // Convenience functions
    tax_t parent(tax_t child) const {
        const auto ind(kh_get(p, pmap_, child));
        return ind == kh_end(pmap_) ? tax_t(-1): kh_val(pmap_, ind);
    }
#define STLFREE(x) do {decltype(x) tmp##x; std::swap(x, tmp##x);} while(0)
    void clear() {
        if(name_map_) khash_destroy(name_map_);
        if(pmap_) khash_destroy(pmap_);
        STLFREE(old_ids_);
        STLFREE(path_map);
        STLFREE(newid_path_map);
    }
#undef STLFREE
    ~TaxonomyReformation() {
        this->clear();
    }
};
using concise_tax_t = TaxonomyReformation;



} // namespace emp

#endif // #ifndef _TX_H__
