#include "test/catch.hpp"
#include "util.h"
#include "circular_buffer.h"

using namespace circ;
using namespace bns;
using namespace std::literals;

// This file tests qmap, encoding, and decoding.


TEST_CASE("circdebuff") {
    SECTION("init") {
        FastCircularQueue<int> buf(16);
        REQUIRE(buf.capacity() == 31);
        buf.emplace_back(137);
        buf.emplace_back(4);
        {
            std::vector<int> vals(&*buf.begin(), &*buf.end());
            REQUIRE(vals.size() == 2);
            REQUIRE(vals[0] == 137);
            REQUIRE(vals[1] == 4);
        }
        buf.emplace_back(13);
        while(buf.size() != buf.capacity()) buf.emplace_back(buf.size());
        buf.pop();
        buf.push_back(13);
        {
            std::vector<int> vals;
            for(const auto val: buf) vals.emplace_back(val);
            REQUIRE(vals.size() == 31);
            REQUIRE(vals[0] == 4);
            REQUIRE(vals[1] == 13);
        }
        buf.pop();
        buf.pop();
        buf.push_back(buf.front());
        REQUIRE(buf.back() == buf.front());
        auto save_front = buf.front();
        REQUIRE(buf.size() == 30);
        buf.push_back(13);
        buf.push_back(14);
        REQUIRE(buf.capacity() == 63);
        REQUIRE(buf.back() == 14);
        REQUIRE(buf.front() == save_front);
    }
}
