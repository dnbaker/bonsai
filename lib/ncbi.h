#ifndef _NBCI_H_
#define _NBCI_H_

#include "htslib/khash.h"
#include "util.h"

namespace kpg {

KHASH_MAP_INIT_INT(p, uint32_t)

khash_t(p) *build_parent_map(const char *fn);
uint32_t lca(khash_t(p) *map, uint32_t a, uint32_t b);

}

#endif // _NBCI_H_
