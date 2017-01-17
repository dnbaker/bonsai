#include "test/catch.hpp"

#include "lib/tx.h"
#include "lib/hashutil.h"
using namespace emp;

TEST_CASE("tax") {
    Taxonomy tax("ref/nodes.dmp", "ref/nameidmap.txt");
    tax.write("test_tax.tx");
    Taxonomy tax2("test_tax.tx");
    //REQUIRE(tax == tax2);
    LOG_DEBUG("Got here\n");
    spvec_t v{3, 7, 1, 0, 4};
    while(v.size() < 30) v.push_back(0);
    std::vector<std::string> paths {
        "ref/viral/GCF_000819615.1_ViralProj14015_genomic.fna.gz",
        "ref/viral/GCF_000820355.1_ViralMultiSegProj14361_genomic.fna.gz",
        "ref/viral/GCF_000820495.2_ViralMultiSegProj14656_genomic.fna.gz",
        "ref/viral/GCF_000836805.1_ViralProj14012_genomic.fna.gz",
        "ref/viral/GCF_000836825.1_ViralProj14017_genomic.fna.gz",
        "ref/viral/GCF_000836845.1_ViralProj14021_genomic.fna.gz",
        "ref/viral/GCF_000836865.1_ViralProj14025_genomic.fna.gz",
        "ref/viral/GCF_000836885.1_ViralProj14030_genomic.fna.gz",
        "ref/viral/GCF_000836905.1_ViralProj14035_genomic.fna.gz",
        "ref/viral/GCF_000836925.1_ViralProj14039_genomic.fna.gz",
        "ref/viral/GCF_000836945.1_ViralProj14044_genomic.fna.gz",
        "ref/viral/GCF_000836965.1_ViralProj14048_genomic.fna.gz",
        "ref/viral/GCF_000836985.1_ViralProj14054_genomic.fna.gz",
        "ref/viral/GCF_000837005.1_ViralProj14058_genomic.fna.gz",
        "ref/viral/GCF_000837025.1_ViralProj14062_genomic.fna.gz",
        "ref/viral/GCF_000837045.1_ViralProj14067_genomic.fna.gz",
        "ref/viral/GCF_000837065.1_ViralProj14071_genomic.fna.gz",
        "ref/viral/GCF_000837085.1_ViralProj14075_genomic.fna.gz",
        "ref/viral/GCF_000837105.1_ViralMultiSegProj14079_genomic.fna.gz",
        "ref/viral/GCF_000837125.1_ViralProj14084_genomic.fna.gz",
        "ref/viral/GCF_000837145.1_ViralProj14089_genomic.fna.gz",
        "ref/viral/GCF_000837165.1_ViralProj14093_genomic.fna.gz",
        "ref/viral/GCF_000837185.1_ViralProj14097_genomic.fna.gz",
        "ref/viral/GCF_000837205.1_ViralMultiSegProj14101_genomic.fna.gz",
        "ref/viral/GCF_000837225.1_ViralProj14105_genomic.fna.gz",
        "ref/viral/GCF_000837245.1_ViralProj14109_genomic.fna.gz",
        "ref/viral/GCF_000837265.1_ViralProj14113_genomic.fna.gz",
        "ref/viral/GCF_000837285.1_ViralMultiSegProj14117_genomic.fna.gz",
        "ref/viral/GCF_000837305.1_ViralMultiSegProj14121_genomic.fna.gz",
        "ref/viral/GCF_000837325.1_ViralProj14125_genomic.fna.gz",
        "ref/viral/GCF_000837345.1_ViralProj14129_genomic.fna.gz",
        "ref/viral/GCF_000837365.1_ViralProj14133_genomic.fna.gz",
        "ref/viral/GCF_000837385.1_ViralProj14137_genomic.fna.gz",
        "ref/viral/GCF_000837405.1_ViralProj14141_genomic.fna.gz",
        "ref/viral/GCF_000837425.1_ViralProj14146_genomic.fna.gz",
        "ref/viral/GCF_000837445.1_ViralProj14150_genomic.fna.gz",
        "ref/viral/GCF_000837465.1_ViralProj14154_genomic.fna.gz",
        "ref/viral/GCF_000837485.1_ViralProj14160_genomic.fna.gz",
        "ref/viral/GCF_000837505.1_ViralMultiSegProj14164_genomic.fna.gz",
        "ref/viral/GCF_000837525.1_ViralProj14168_genomic.fna.gz"
    };
    LOG_DEBUG("Got here\n");
    Spacer sp(31, 31, v);
    LOG_DEBUG("Making set\n");
    kgset_t set(paths, sp);
    REQUIRE(set.size() == paths.size());
    count::Counter<std::vector<std::uint64_t>> counts(bitmap_t(set).to_counter());
    LOG_DEBUG("Weight: %zu. Number of bit patterns: %zu\n", set.weight(), counts.size());
    for(auto &i: counts) {
        LOG_DEBUG("Number of 64-bit integers used: %zu\n", i.first.size());
        for(auto j: i.first) {
            for(auto k(0); k < CHAR_BIT * sizeof(j); ++k) {
                fputc('0' + (j & 1), stderr);
                j >>= 1;
            }
        }
        fputc('\n', stderr);
    }
    constexpr char size_arr[] {"ACADasdadasdaj2387asdfadsfakjsfaksdjhfakjsdhfjkasdf"};
    REQUIRE(spop(size_arr) == vec_popcnt(size_arr, sizeof(size_arr) - 1));
    LOG_DEBUG("Value of thing is");
    const uint64_t SEED[4] = {0x818c3f78ull, 0x672f4a3aull, 0xabd04d69ull,
                              0x12b51f95ull};
    const unsigned __int128 m =
        *reinterpret_cast<const unsigned __int128 *>(&SEED[0]);
    const unsigned __int128 a =
        *reinterpret_cast<const unsigned __int128 *>(&SEED[2]);
    cuckoofilter::putu128(m, stdout);
    cuckoofilter::putu128(a, stdout);
}
