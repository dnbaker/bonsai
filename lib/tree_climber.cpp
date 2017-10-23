#include "lib/tree_climber.h"
#include <stdexcept>
#include <system_error>

namespace emp { namespace tree {


std::string make_fname(const Database<khash_t(c)> &db, const std::string &fld, const tax_t tax) {
    char buf[128];
    std::sprintf(buf, "%s%u.w%u.k%u.kmers.bin", fld.data(), tax, db.w_, db.k_);
    return buf;
}

std::pair<std::vector<std::string>, std::unordered_set<tax_t>> invert_lca_map(const Database<khash_t(c)> &db, const char *folder, int prebuilt) {
    std::vector<std::string> ret;
    std::unordered_set<tax_t> parents;
    std::unordered_map<tax_t, std::FILE *> ofps;
    std::unordered_map<tax_t, std::uint64_t> ofp_counts;
    khash_t(c) *map(db.db_); // Does not own; shorthand.
    assert(map);
    std::string fld;
    int fwritten;
#if !NDEBUG
    std::size_t n_processed(0);
#endif
    if(folder && *folder) fld = folder, fld += '/';
    if(prebuilt == 0) {
        for(khiter_t ki(0); ki != kh_end(map); ++ki) {
            if(!kh_exist(map, ki)) continue;
            auto m(ofps.find(kh_val(map, ki)));
            auto mc(ofp_counts.find(kh_val(map, ki)));
            if(m == ofps.end()) {
                LOG_DEBUG("Adding %u to map\n", kh_val(map, ki));
                ret.emplace_back(make_fname(db, fld, kh_val(map, ki)));
                m = ofps.emplace(kh_val(map, ki), std::fopen(ret[ret.size() - 1].data(), "wb")).first;
                if(m->second == nullptr) {
                    char buf2[512];
                    std::sprintf(buf2, "Could not open file at path %s\n", ret[ret.size() - 1].data());
                    perror(buf2);
                    exit(1);
                }
                LOG_DEBUG("Opened file at %s\n", ret[ret.size() - 1].data());
                mc = ofp_counts.emplace(kh_val(map, ki), 1).first;
                parents.insert(kh_val(map, ki));
                std::uint64_t tmp(-1);
                fwritten = std::fwrite(&tmp, sizeof(tmp), 1, m->second);
                if(fwritten != 1) {
                    char buf2[256];
                    std::sprintf(buf2, "Could not write dummy value. Size written: %i\n", fwritten);
                    perror(buf2), exit(1);
                }
            } else ++mc->second;
            if((fwritten = std::fwrite(&kh_key(map, ki), sizeof(kh_key(map, ki)), 1, m->second)) != 1) {
                char buf2[256];
                std::sprintf(buf2, "Could not write dummy value. Size written: %i. taxid: %u\n", fwritten, kh_val(map, ki));
                perror(buf2), exit(1);
            }
#if !NDEBUG
            if((++n_processed & 0xFFFFFF) == 0) LOG_DEBUG("Processed %zu so far\n", n_processed);
#endif
        }
        for(auto &pair: ofps) {
            std::rewind(pair.second);
            std::fwrite(&ofp_counts[pair.first], 1, sizeof(std::uint64_t), pair.second);
            std::fclose(pair.second);
        }
    } else {
        LOG_DEBUG("Just getting filenames\n");
        std::unordered_set<tax_t> taxes;
        taxes.reserve(1 << 12);
        for(khiter_t ki(0); ki != kh_end(map); ++ki) {
            if(!kh_exist(map, ki)) continue;
            auto m(taxes.find(kh_val(map, ki)));
            if(m == taxes.end()) {
                ret.emplace_back(make_fname(db, fld, kh_val(map, ki)));
                taxes.insert(kh_val(map, ki));
            }
        }
    }
    return std::pair<std::vector<std::string>, std::unordered_set<tax_t>>(std::move(ret), std::move(parents));
}

khash_t(p) *pruned_taxmap(std::vector<std::string> &paths, khash_t(p) *taxmap, khash_t(name) *name_hash) {
    std::set<tax_t> found, keep;
    std::uint64_t ki1, ki2;
    for(auto &path: paths) found.insert(get_taxid(path.data(), name_hash));
    {
        auto m(found.find(-1u));
        if(m != found.end()) found.erase(m);
    }
    for(auto i: found) {
        keep.insert(i);
        ki1 = kh_get(p, taxmap, i);
        if(ki1 == kh_end(taxmap)) {
            LOG_DEBUG("tax %u missing from map. Skipping in tree pruning.\n", i);
            continue;
        }
        while((i = kh_val(taxmap, ki1))) {
            keep.insert(i);
            if((ki1 = kh_get(p, taxmap, i)) == kh_end(taxmap)) {
                LOG_EXIT("Missing taxid for %u. ZOMGZZZZ Skipping....", i);
                break;
            }
        }
    }
    int khr;
    khash_t(p) *ret(kh_init(p));
    for(const auto i: keep) {
        ki1 = kh_put(p, ret, i, &khr);
        ki2 = kh_get(p, taxmap, i);
        if(ki2 == kh_end(taxmap)) continue;
        kh_val(ret, ki1) = kh_val(taxmap, ki2);
    }
    return ret;
}

khash_t(all) *load_binary_kmerset(const char *path) {
    std::FILE *fp(fopen(path, "rb"));
    if(fp == nullptr) throw std::system_error(std::error_code(2, std::system_category()), std::string("Cannot open path at ") + path + ".\n");
    khash_t(all) *ret(kh_init(all));
    std::uint64_t n;
    std::fread(&n, 1, sizeof(n), fp);
    if(kh_resize(all, ret, n) < 0) LOG_EXIT("Could not resize hash table to next power of 2 above %zu. New size: %zu\n", n, kh_n_buckets(ret));
    LOG_DEBUG("About to place %zu elements into a hash table of max size %zu\n", n, kh_n_buckets(ret));
    for(int khr; std::fread(&n, 1, sizeof(std::uint64_t), fp) == sizeof(std::uint64_t); kh_put(all, ret, n, &khr));
    std::fclose(fp);
#if !NDEBUG
    for(khiter_t ki(0); ki < kh_end(ret); ++ki) {
        if(kh_exist(ret, ki)) assert(kh_get(all, ret, kh_key(ret, ki)) != kh_end(ret));
    }
    fp = fopen(path, "rb");
    std::fread(&n, 1, sizeof(std::uint64_t), fp); // Skip first number.
    while(std::fread(&n, 1, sizeof(std::uint64_t), fp) == sizeof(std::uint64_t)) assert(kh_get(all, ret, n) != kh_end(ret));
    std::fclose(fp);
#endif
    return ret;
}


bitvec_t load_binary_kmers(const char *path) {
    std::FILE *fp(fopen(path, "rb"));
    bitvec_t ret;
    std::uint64_t n, ind(0);
    std::fread(&n, sizeof(n), 1, fp);
    ret.resize(n, false); // Initializes to 0 unnecessarily. Better than passing it to a temporary variable every time, I think.
    while(std::fread(ret.data() + ind++, sizeof(std::uint64_t), 1, fp) == sizeof(std::uint64_t));
    std::fclose(fp);
    return ret;
}

} /* namespace tree */ } // namespace emp
