// RollingHasherSet must give hasher i the seeds 2i and 2i + 1 drawn from
// std::mt19937_64(seedseed), in that order, on every platform, so that the same
// seedseed gives the same hashes everywhere.
#include "test/catch.hpp"
#include "encoder.h"
#include <random>
#include <vector>
using namespace bns;
TEST_CASE("RollingHasherSet draws each hasher's seeds in order") {
    const std::vector<unsigned> kvals{21, 33, 40};
    RollingHasherSet<uint64_t> rhs(kvals, true, DNA, 1337u);
    std::mt19937_64 mt(1337u);
    for(size_t hi = 0; hi < kvals.size(); ++hi) {
        const uint64_t seed1 = mt();
        const uint64_t seed2 = mt();
        REQUIRE(rhs.hashers_[hi].seed1_ == seed1);
        REQUIRE(rhs.hashers_[hi].seed2_ == seed2);
    }
}
