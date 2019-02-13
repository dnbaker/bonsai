#include "tx.h"
#include "test/catch.hpp"

using namespace bns;

std::vector<std::string> PATHS {
    "ec/GCF_000005845.2_ASM584v2_genomic.fna.gz",
    "ec/GCF_000007445.1_ASM744v1_genomic.fna.gz",
    "ec/GCF_000008865.1_ASM886v1_genomic.fna.gz",
    "ec/GCF_000009565.1_ASM956v1_genomic.fna.gz",
    "ec/GCF_000010245.2_ASM1024v1_genomic.fna.gz",
    "ec/GCF_000010385.1_ASM1038v1_genomic.fna.gz",
    "ec/GCF_000010485.1_ASM1048v1_genomic.fna.gz",
    "ec/GCF_000010745.1_ASM1074v1_genomic.fna.gz",
    "ec/GCF_000010765.1_ASM1076v1_genomic.fna.gz",
    "ec/GCF_000013265.1_ASM1326v1_genomic.fna.gz",
    "ec/GCF_000013305.1_ASM1330v1_genomic.fna.gz",
    "ec/GCF_000014845.1_ASM1484v1_genomic.fna.gz",
    "ec/GCF_000017745.1_ASM1774v1_genomic.fna.gz",
    "ec/GCF_000017765.1_ASM1776v1_genomic.fna.gz",
    "ec/GCF_000017985.1_ASM1798v1_genomic.fna.gz",
    "ec/GCF_000019385.1_ASM1938v1_genomic.fna.gz",
    "ec/GCF_000019425.1_ASM1942v1_genomic.fna.gz",
    "ec/GCF_000019645.1_ASM1964v1_genomic.fna.gz",
    "ec/GCF_000021125.1_ASM2112v1_genomic.fna.gz",
    "ec/GCF_000022225.1_ASM2222v1_genomic.fna.gz",
    "ec/GCF_000022345.1_ASM2234v1_genomic.fna.gz",
    "ec/GCF_000022665.1_ASM2266v1_genomic.fna.gz",
    "ec/GCF_000465235.1_ASM46523v1_genomic.fna.gz",
    "ec/GCF_000494775.1_ASM49477v1_genomic.fna.gz",
    "ec/GCF_000583755.1_ASM58375v1_genomic.fna.gz",
    "ec/GCF_000583775.1_ASM58377v1_genomic.fna.gz",
    "ec/GCF_000583795.1_ASM58379v1_genomic.fna.gz",
    "ec/GCF_000954195.1_ASM95419v1_genomic.fna.gz",
    "ec/GCF_001305715.1_ASM130571v1_genomic.fna.gz",
    "ec/GCF_001417635.1_ASM141763v1_genomic.fna.gz",
    "ec/GCF_001483845.1_ASM148384v1_genomic.fna.gz",
    "ec/GCF_001639125.1_ASM163912v1_genomic.fna.gz",
    "ec/GCF_001717605.1_ASM171760v1_genomic.fna.gz",
    "ec/GCF_001865455.1_ASM186545v1_genomic.fna.gz",
    "ec/GCF_001865475.1_ASM186547v1_genomic.fna.gz",
    "ec/GCF_001865495.1_ASM186549v1_genomic.fna.gz",
    "ec/GCF_001865515.1_ASM186551v1_genomic.fna.gz",
    "ec/GCF_001865535.1_ASM186553v1_genomic.fna.gz",
    "ec/GCF_001865555.1_ASM186555v1_genomic.fna.gz",
    "ec/GCF_001936355.1_ASM193635v1_genomic.fna.gz"
};

const char *TAX_PATH  = "../ref/nodes.dmp";
const char *NAME_PATH = "../ref/nameidmap.txt";

TEST_CASE("tx") {
    khash_t(p) *old_tax(build_parent_map(TAX_PATH));
    khash_destroy(old_tax);
}
