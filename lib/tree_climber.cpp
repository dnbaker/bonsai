#include "lib/tree_climber.h"

namespace emp { namespace tree {

std::vector<std::string> invert_lca_map(Database<khash_t(c)> &db, const char *folder, int prebuilt) {
    std::vector<std::string> ret;
    std::unordered_map<std::uint32_t, std::FILE *> ofps;
    std::unordered_map<std::uint32_t, std::uint64_t> ofp_counts;
    khash_t(c) *map(db.db_); // Does not own; shorthand.
    assert(map);
    std::string fld;
    int fwritten;
    if(folder && *folder) fld = folder, fld += '/';
    if(prebuilt == 0) {
        for(khiter_t ki(0); ki != kh_end(map); ++ki) {
            if(!kh_exist(map, ki)) continue;
            auto m(ofps.find(kh_val(map, ki)));
            auto mc(ofp_counts.find(kh_val(map, ki)));
            if(m == ofps.end()) {
                LOG_DEBUG("Adding %u to map\n", kh_val(map, ki));
                const std::uint64_t tmp(UINT64_C(-1));
                char buf[256];
                std::sprintf(buf, "%s%u.kmers.bin", fld.data(), kh_val(map, ki));
                ret.emplace_back(buf);
                m = ofps.emplace(kh_val(map, ki), std::fopen(buf, "wb")).first;
                if(m->second == nullptr) {
                    char buf2[512];
                    std::sprintf(buf2, "Could not open file at path %s\n", buf);
                    perror(buf2);
                    exit(1);
                }
                LOG_DEBUG("Opened file at %s\n", buf);
                mc = ofp_counts.emplace(kh_val(map, ki), 1).first;
                fwritten = std::fwrite(&tmp, sizeof(tmp), 1, m->second);
#if 0
                if(fwritten != sizeof(tmp)) {
                    char buf2[256];
                    std::sprintf(buf2, "Could not write dummy value. Size written: %i\n", fwritten);
                    perror(buf), exit(1);
                }
                LOG_DEBUG("Wrote one entry!\n");
#endif
            } else ++mc->second;
            fwritten = std::fwrite(&kh_key(map, ki), sizeof(kh_key(map, ki)), 1, m->second);
#if 0
            if(fwritten != sizeof(kh_key(map, ki))) {
                char buf2[256];
                std::sprintf(buf2, "Could not write dummy value. Size written: %i. taxid: %u\n", fwritten, kh_val(map, ki));
                perror(buf2), exit(1);
            }
#endif
#if !NDEBUG
            if((ki & 0xFFFFFF) == 0) LOG_DEBUG("Processed %zu so far\n", (std::size_t)ki);
#endif
        }
        for(auto &pair: ofps) {
            std::rewind(pair.second);
            std::fwrite(&ofp_counts[pair.first], 1, sizeof(std::uint64_t), pair.second);
            std::fclose(pair.second);
        }
    } else {
        for(khiter_t ki(0); ki != kh_end(map); ++ki) {
            if(!kh_exist(map, ki)) continue;
            auto m(ofps.find(kh_val(map, ki)));
            if(m == ofps.end()) {
                char buf[256];
                std::sprintf(buf, "%s%u.kmers.bin", fld.data(), kh_val(map, ki));
                ret.emplace_back(buf);
            }
        }
    }
    return ret;
}

khash_t(p) *pruned_taxmap(std::vector<std::string> &paths, khash_t(p) *taxmap, khash_t(name) *name_hash) {
    std::set<std::uint32_t> found, keep;
    std::uint64_t ki1, ki2;
    for(auto &path: paths) found.insert(get_taxid(path.data(), name_hash));
    for(auto i: found) {
        keep.insert(i);
        ki1 = kh_get(p, taxmap, i);
        i = kh_val(taxmap, ki1);
        while(i) {
            keep.insert(i);
            ki1 = kh_get(p, taxmap, i);
            i = kh_val(taxmap, ki1);
        }
    }
    int khr;
    khash_t(p) *ret(kh_init(p));
    for(const auto i: keep) {
        ki1 = kh_put(p, ret, i, &khr);
        ki2 = kh_get(p, taxmap, i);
        kh_val(ret, ki1) = kh_val(taxmap, ki2);
    }
    return ret;
}

khash_t(all) *load_binary_kmerset(const char *path) {
    std::FILE *fp(fopen(path, "rb"));
    khash_t(all) *ret(kh_init(all));
    std::uint64_t n;
    int khr;
    std::fread(&n, 1, sizeof(n), fp);
    kh_resize(all, ret, n);
    while(std::fread(&n, 1, sizeof(std::uint64_t), fp) == sizeof(std::uint64_t))
        kh_put(all, ret, n, &khr);
    std::fclose(fp);
    return ret;
}


std::vector<std::uint64_t> load_binary_kmers(const char *path) {
    std::FILE *fp(fopen(path, "rb"));
    std::vector<std::uint64_t> ret;
    std::uint64_t n, ind(0);
    std::fread(&n, sizeof(n), 1, fp);
    ret.resize(n); // Initializes to 0 unnecessarily. Better than passing it to a temporary variable every time, I think.
    while(std::fread(ret.data() + ind++, sizeof(std::uint64_t), 1, fp) == sizeof(std::uint64_t));
    std::fclose(fp);
    return ret;
}

} /* namespace tree */ } // namespace emp
