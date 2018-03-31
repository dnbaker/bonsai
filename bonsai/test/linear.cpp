#include "test/catch.hpp"

#include <set>
#include "tx.h"
#include "bitmap.h"
#include "linear/linear.h"
using namespace bns;

TEST_CASE("set") {
    linear::set<uint32_t> linset;
    for(size_t i(0); i < 20; linset.insert(i++));
    std::set<uint32_t> oset;
    for(const auto el: linset) {
        oset.insert(el);
        REQUIRE(linset.contains(el));
    }
    for(size_t i(0); i < 20; ++i) {
        REQUIRE(oset.find(i) != oset.end());
        REQUIRE(linset.find(i) != linset.end());
    }
    assert(linset.size() == oset.size());
}

TEST_CASE("counter") {
    linear::counter<uint32_t> linset;
    for(size_t i(0); i < 20; ++i) for(size_t j(0); j < i; ++j, linset.add(i));
    for(size_t i(0); i < 20; ++i) REQUIRE(linset.count(i) == i);
}
