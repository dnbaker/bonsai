#include "tree_query.h"

namespace emp {
#if 0
KrakenLcaMap::KrakenLcaMap(const char *dbpath, const char *idxpath) {
    kraken::QuickFile qf;
    qf.open_file(dbpath);
    qf.load_file(dbpath);
}

KrakenLcaMap::~KrakenLcaMap() {

}
#endif

KmerLcaMap::~KmerLcaMap() {
}

u32 EmpLcaMap::get_lca(u64 kmer) {
    return db_.get_lca(kmer);
}

EmpLcaMap::~EmpLcaMap() {
}

} // namespace emp
