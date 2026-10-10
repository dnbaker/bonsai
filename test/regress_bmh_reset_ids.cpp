#include "test/catch.hpp"
#include "sketch/bmh.h"
#include <cstdint>

// reset() must clear the sampled ids and counts of a BagMinHash that tracks
// them, so a sketch reused for an empty input does not report the ids and
// counts of the previous input.

TEST_CASE("BagMinHash reset clears the sampled ids and counts") {
    const size_t m = 16;
    sketch::BagMinHash2<double> h(m, true, true);
    for(uint64_t id = 1; id <= 500; ++id) h.update(id, double(1 + id % 3));
    h.finalize();
    size_t nset = 0;
    for(size_t i = 0; i < m; ++i) nset += h.ids()[i] != 0;
    REQUIRE(nset == m);
    h.reset();
    h.finalize();
    size_t nstale = 0;
    for(size_t i = 0; i < m; ++i) nstale += h.ids()[i] != 0 || h.idcounts()[i] != 0.;
    REQUIRE(nstale == 0);
}
