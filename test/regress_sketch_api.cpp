// Sketch merge API symmetry: BBitMinHasher provides operator+ matching
// operator+= (as hllbase_t does), and the copy-merge equals in-place merge.
#include "test/catch.hpp"

#include "sketch/bbmh.h"
#include "util.h"

using namespace bns;
using namespace sketch;

TEST_CASE("BBitMinHasher operator+ merges like operator+=") {
    BBitMinHasher<u64> a(11, 16), b(11, 16);
    for(u64 i(0); i < 5000; ++i) {
        a.addh(i * 2);
        b.addh(i * 2 + 1);
    }
    BBitMinHasher<u64> merged(a + b);
    BBitMinHasher<u64> in_place(a);
    in_place += b;
    auto fm(merged.finalize()), fi(in_place.finalize());
    REQUIRE(fm == fi);
    REQUIRE(fm.jaccard_index(fi) == 1.0);
}
