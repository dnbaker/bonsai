#include "catch.hpp"
#include "hll/hll.h"
#include "bonsai/include/util.h"

template<typename T> class TD;

TEST_CASE("SIMDhash") {
    hll::WangHash hf;
    hll::MurFinHash mfh;
    std::mt19937_64 mt;
    __attribute__((aligned(alignof(vec::SIMDTypes<uint64_t>::ValueType)))) u64 arr[vec::SIMDTypes<uint64_t>::COUNT];
    const size_t ntrials = 100000;
    __attribute__((aligned(alignof(vec::SIMDTypes<uint64_t>::ValueType)))) u64 arrret[vec::SIMDTypes<uint64_t>::COUNT];
    for(size_t i(0); i < ntrials; ++i) {
        u64 val = mt();
        std::fill(std::begin(arr), std::end(arr), val);
        //std::fill(std::begin(arr), std::end(arr), val);
        auto hv64 = hf(arr[0]);
        auto hv256 = hf(*(vec::SIMDTypes<uint64_t>::Type *)arr);
        std::memcpy(arrret, &hv256, sizeof(hv256));
        REQUIRE(hv64 == arrret[0]);
        REQUIRE(hv64 == arrret[1]);
#if __AVX2__ || HAS_AVX_512
        REQUIRE(hv64 == arrret[2]);
        REQUIRE(hv64 == arrret[3]);
#endif
        std::fill(std::begin(arr), std::end(arr), val);
        hv64 = mfh(arr[0]);
        hv256 = mfh(*(vec::SIMDTypes<uint64_t>::Type *)arr);
        std::memcpy(arrret, &hv256, sizeof(hv256));
        REQUIRE(hv64 == arrret[0]);
        REQUIRE(hv64 == arrret[1]);
#if __AVX2__ || HAS_AVX_512
        REQUIRE(hv64 == arrret[2]);
        REQUIRE(hv64 == arrret[3]);
#endif
    }
}
