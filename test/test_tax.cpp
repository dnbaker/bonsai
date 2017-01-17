#include "test/catch.hpp"

#include "lib/tx.h"
#include "lib/hashutil.h"
using namespace emp;

TEST_CASE("tax") {
    Taxonomy tax("ref/nodes.dmp", "ref/nameidmap.txt");
    tax.write("test_tax.tx");
    Taxonomy tax2("test_tax.tx");
    REQUIRE(tax == tax2);
    LOG_DEBUG("Got here\n");
    spvec_t v{3, 7, 1, 0, 4};
    while(v.size() < 30) v.push_back(0);
    std::vector<std::string> paths;
    {
        kstring_t ks{0, 0, 0};
        const char cmd[] {"ls ref/viral/ | grep fna.gz | head -n 500"};
        fprintf(stderr, "cmd: %s\n", cmd);
        FILE *fp(popen(cmd, "r"));
        char buf[512];
        while(fgets(buf, sizeof buf, fp)) {
            std::string tmp("ref/viral/");
            tmp += buf;
            tmp.pop_back();
            paths.emplace_back(std::move(tmp));
        }
        fclose(fp);
    }
    LOG_DEBUG("Got here\n");
    Spacer sp(31, 31, v);
    LOG_DEBUG("Making set\n");
    kgset_t set(paths, sp);
    REQUIRE(set.size() == paths.size());
    count::Counter<std::vector<std::uint64_t>> counts(bitmap_t(set).to_counter());
    LOG_DEBUG("Weight: %zu. Number of bit patterns: %zu. Total weight of all bit patterns: %zu\n", set.weight(), counts.size(), counts.total());
    counts.print_counts(stderr);
#if 0
    for(auto &i: counts) {
        for(auto j: i.first) {
            for(auto k(0u); k < CHAR_BIT * sizeof(j); ++k) {
                fputc('0' + (j & 1), stderr);
                j >>= 1;
            }
        }
        fputc('\n', stderr);
    }
#endif
    constexpr char size_arr[] {"ACADasdadasdaj2387asdfadsfakjsfaksdjhfakjsdhfjkasdasdfsdasddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddafasdfasdfasdf"};
    REQUIRE(spop(size_arr) == vec_popcnt(size_arr, sizeof(size_arr)));
#if 0
    unsigned __int128 i(((unsigned __int128)0x12b51f95ull << 64) | 0xabd04d69ull),
                      j(((unsigned __int128)0x672f4a3aull << 64) | 0x818c3f78ull);
    cuckoofilter::putu128(i, stdout);
    cuckoofilter::putu128(j, stdout);
#endif
}
