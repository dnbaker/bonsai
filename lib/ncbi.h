#ifndef _NBCI_H_
#define _NBCI_H_

#include "htslib/khash.h"
#include "util.h"

namespace kpg {

khash_t(p) *build_parent_map(const char *fn);

}

#endif // _NBCI_H_
