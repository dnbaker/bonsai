#ifndef __DB_H__
#define __DB_H__

#include "feature_min.h"
#include "cuckoohash_map.hh"
#include "city_hasher.hh"

namespace kpg {
typedef cuckoohash_map<uint64_t, uint32_t, CityHasher<uint64_t>> chm_t;

template<uint64_t (*score)(uint64_t, void *)>
void build_minimized_database(khash_t(64) *td_map, const Spacer &sp, std::vector<std::string> &paths, chm_t &ret);
extern std::map<int, void (*)(khash_t(64) *, const Spacer &, std::vector<std::string> &, chm_t &)> dbbuild_map;

} // namespace kpg

#endif // #ifndef db
