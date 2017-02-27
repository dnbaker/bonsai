#include "test/catch.hpp"

#include "lib/tx.h"
#include "lib/bitmap.h"
#include "lib/tree_climber.h"
using namespace emp;

TEST_CASE("tax") {
    Taxonomy tax("ref/nodes.dmp", "ref/nameidmap.txt");
    tax.write("test_tax.tx");
    Taxonomy tax2("test_tax.tx");
    REQUIRE(tax == tax2);
    spvec_t v{3, 7, 1, 0, 4};
    while(v.size() < 12) v.push_back(0);
    std::vector<std::string> paths;
    {
        const char cmd[] {"ls ref/viral/ | grep fna.gz | head -n 40"};
        std::FILE *fp(popen(cmd, "r"));
        if(fp == nullptr) throw "a party";//throw std::system_error(1, std::string("Could not call command '") + cmd + "'\n");
        int c;
        size_t ind(0);
        ssize_t len;
        size_t size;
        char *buf(nullptr);
        const std::string prefix("ref/viral/");
        while((len = getline(&buf, &size, fp)) >= 0) {
            buf[len - 1] = '\0';
            paths.emplace_back(prefix + buf);
        }
        free(buf);
        std::fclose(fp);
    }
    Spacer sp(13, 13, v);
    kgset_t set(paths, sp);
    REQUIRE(set.size() == paths.size());
    count::Counter<std::vector<std::uint64_t>> counts(bitmap_t(set).to_counter());
    adjmap_t adj(counts);
    //LOG_DEBUG("Weight: %zu. Number of bit patterns: %zu. Total weight of all bit patterns: %zu\n", set.weight(), counts.size(), counts.total());

    //counts.print_counts(stderr);
    //counts.print_hist(stderr);
    std::vector<std::uint64_t> thing;
    for(size_t i(0); i < 1 << 16; ++i) thing.push_back(((uint64_t)rand() << 32) | rand());
    REQUIRE(spop(thing) == popcnt::vec_popcnt(thing));
}

TEST_CASE("bitstrings") {
    std::vector<std::uint64_t> v1, v2;
    for(auto i: {1,177,123232,1222, 3344411, 11232}) v1.emplace_back(i), v2.emplace_back(i);
    count::Counter<std::vector<std::uint64_t>> counter;
    counter.add(v1);
    counter.add(v2);
    counter.print_vec();
}
