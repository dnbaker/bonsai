#ifndef _NBCI_H_
#define _NBCI_H_

#include "htslib/khash.h"
#include "util.h"

namespace kpg {

KHASH_MAP_INIT_INT(p, uint32_t)

khash_t(p) *build_parent_map(const char *fn);

}

#endif // _NBCI_H_
