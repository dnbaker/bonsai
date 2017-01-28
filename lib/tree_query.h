#ifndef TREE_QUERY_H__
#define TREE_QUERY_H__
#include "kraken/src/krakenutil.hpp"
#include "kraken/src/krakendb.hpp"
#include "kmerutil.h"
#include "tx.h"
#include "database.h"

namespace emp {

class KmerLcaMap {
public:
    virtual std::uint32_t get_lca(std::uint64_t) = 0;
    virtual ~KmerLcaMap() = 0;
};

#if 0
class KrakenLcaMap: KmerLcaMap {
    kraken::KrakenDB      db_;
    kraken::KrakenDBIndex dbi_;
    kraken::KmerScanner   ks_;
public:
    virtual std::uint32_t get_lca(std::uint64_t) override;
    virtual ~KrakenLcaMap();
};
#endif

class EmpLcaMap: KmerLcaMap {
    //Taxonomy tx_;
    Database<khash_t(c)> db_;
public:
    EmpLcaMap(): db_("LoadFromThisNotRealFile") {
    }
    virtual std::uint32_t get_lca(std::uint64_t) override;
    virtual ~EmpLcaMap();
};

}

#endif // #ifndef TREE_QUERY_H__
