#include "tx.h"
#include "test/catch.hpp"

using namespace bns;

std::vector<std::string> PATHS {
    "ec/GCF_000005845.2_ASM584v2_genomic.fna.gz",
    "ec/GCF_000007445.1_ASM744v1_genomic.fna.gz",
    "ec/GCF_000008865.1_ASM886v1_genomic.fna.gz",
    "ec/GCF_000009565.1_ASM956v1_genomic.fna.gz",
    "ec/GCF_000010245.2_ASM1024v1_genomic.fna.gz",
};

const char *TAX_PATH  = "../ref/nodes.dmp";
const char *NAME_PATH = "../ref/nameidmap.txt";

TEST_CASE("tx") {
    if(!bns::isfile(TAX_PATH)) {
        std::fprintf(stderr, "Tax path %s does not exist. Skipping test.\n", TAX_PATH);
    } else {
        khash_t(p) *old_tax(build_parent_map(TAX_PATH));
        khash_destroy(old_tax);
    }
}
