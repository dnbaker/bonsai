#include "lib/tree_climber.h"

namespace emp { namespace tree {

std::vector<std::string> invert_lca_map(Database<khash_t(c)> &db, const char *folder) {
    std::vector<std::string> ret;
    std::unordered_map<std::uint32_t, std::FILE *> ofps;
    std::unordered_map<std::uint32_t, std::uint64_t> ofp_counts;
    khash_t(c) *map(db.db_); // Does not own; shorthand.
    for(khiter_t ki(0); ki != kh_end(map); ++ki) {
        if(!kh_exist(map, ki)) continue;
        auto m(ofps.find(kh_val(map, ki)));
        auto mc(ofp_counts.find(kh_val(map, ki)));
        if(m == ofps.end()) {
            char buf[64];
            const std::uint64_t tmp(UINT64_C(-1));
            std::sprintf(buf, "%s/%u.kmers.bin", folder, kh_val(map, ki));
            ret.emplace_back(buf);
            m = ofps.emplace(kh_val(map, ki), std::fopen(buf, "wb")).first;
            mc = ofp_counts.emplace(kh_val(map, ki), 0).first;
            std::fwrite(&tmp, 1, sizeof(tmp), m->second);
        }
        std::fwrite(&kh_key(map, ki), sizeof(kh_key(map, ki)), 1, m->second);
        ++mc->second;
    }
    for(auto &pair: ofps) {
        std::rewind(pair.second);
        std::fwrite(&ofp_counts[pair.first], 1, sizeof(std::uint64_t), pair.second);
        std::fclose(pair.second);
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

std::vector<std::uint64_t> load_binary_kmers(const char *path) {
    std::FILE *fp(fopen(path, "rb"));
    std::vector<std::uint64_t> ret;
    std::uint64_t n, ind(0);
    std::fread(&n, 1, sizeof(n), fp);
    ret.resize(n); // Initializes to 0 unnecessarily. Better than passing it to a temporary variable every time, I think.
    while(std::fread(ret.data() + ind++, 1, sizeof(std::uint64_t), fp) == sizeof(std::uint64_t));
    std::fclose(fp);
    return ret;
}

} /* namespace tree */ } // namespace emp
