#ifndef _NBCI_H_
#define _NBCI_H_

#include "htslib/khash.h"
#include "util.h"

namespace kpg {

KHASH_MAP_INIT_INT(p, uint32_t)

khash_t(p) *build_parent_map(const char *fn);
void write_parent_map(khash_t(p) *hash, const char *fn);
khash_t(p) *load_parent_map(khash_t(p) *hash, const char *fn);

}

#endif // _NBCI_H_
