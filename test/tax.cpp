#include "test/catch.hpp"

#include "lib/tx.h"
#include "lib/bitmap.h"
using namespace emp;

TEST_CASE("tax") {
    spvec_t v{3, 7, 1, 0, 4};
    while(v.size() < 12) v.push_back(0);
    std::vector<std::string> paths;
    {
        const char cmd[] {"ls ec | grep fna.gz | head -n 40"};
        std::FILE *fp(popen(cmd, "r"));
        if(fp == nullptr) throw "a party";//throw std::system_error(1, std::string("Could not call command '") + cmd + "'\n");
        int c;
        size_t ind(0);
        ssize_t len;
        size_t size;
        char *buf(nullptr);
        const std::string prefix("ec/");
        while((len = getline(&buf, &size, fp)) >= 0) {
            buf[len - 1] = '\0';
            paths.emplace_back(prefix + buf);
            if(!std::ifstream(paths.back()).good()) throw std::runtime_error(ks::sprintf("Could not open path at %s\n", paths.back().data()).data());
        }
        std::free(buf);
        std::fclose(fp);
    }
    std::fprintf(stderr, "Parsed paths. Number of paths: %zu\n", paths.size());
    Spacer sp(13, 13, v);
    kgset_t set(paths, sp);
    REQUIRE(set.size() == paths.size());
    count::Counter<bitvec_t> counts(bitmap_t(set).to_counter());
    adjmap_t adj(counts);
    std::fprintf(stderr, "Made adjmap");

    counts.print_counts(stderr);
    counts.print_hist(stderr);
    bitvec_t thing;
    thing.reserve(1 << 16);
    for(size_t i(0); i < 1 << 16; ++i)
        thing.emplace_back(((uint64_t)rand() << 32) | rand());
    std::fprintf(stderr, "Made thing");
}

TEST_CASE("bitstrings") {
    bitvec_t v1, v2;
    for(auto i: {1,177,123232,1222, 3344411, 11232}) v1.emplace_back(i), v2.emplace_back(i);
    count::Counter<bitvec_t> counter;
    counter.add(v1);
    counter.add(v2);
    counter.print_vec();
}
