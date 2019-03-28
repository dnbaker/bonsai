#include "test/catch.hpp"
#include "spacer.h"
#include "encoder.h"
#include <algorithm>
#include <numeric>

using namespace bns;
using EncType = Encoder<score::Lex>;
using namespace std::literals;

// This file tests qmap, encoding, and decoding.


TEST_CASE( "Spacer encodes and decodes contiguous, unminimized seeds correctly.", "[contiguous]") {
    spvec_t v;
    const char *test_str("ACATGCTAGCATGCTGACTGACTGATCGATCGTA");
    SECTION("Add spaces") {
        std::string test(test_str);
        v = {1, 2}; // Skip 1 base after the first nucleotide, then skip two, then don't skip any.
        test.resize(34); // Just long enough to encode a k-31 with 3 total spaces in it.
        while(v.size() < 30) v.push_back(0);
        Spacer sp(31, 31, v);
        EncType enc(&test[0], 34, sp, nullptr, true);
        REQUIRE(enc.canonicalize());
        REQUIRE(sp.s_.size() == 30);
        REQUIRE(std::accumulate(sp.s_.begin(), sp.s_.end(), 0u) == sp.s_.size() + std::accumulate(v.begin(), v.end(), 0));
        std::string to_use = test;
        to_use[1] = '-'; to_use[3] = '-'; to_use[4] = '-';
        REQUIRE(to_use != test);
        REQUIRE(sp.to_string(enc.next_kmer()) == to_use);
    }
    SECTION("next_minimizer rather than kmer") {
        std::string test(test_str);
        v = {1, 2};
        test.resize(34);
        while(v.size() < 30) v.push_back(0);
        Spacer sp(31, 31, v);
        EncType enc(&test[0], 34, sp, nullptr, true);
        REQUIRE(enc.canonicalize());
        REQUIRE(sp.s_.size() == 30);
        REQUIRE(std::accumulate(sp.s_.begin(), sp.s_.end(), 0u) == sp.s_.size() + std::accumulate(v.begin(), v.end(), 0));
        std::string to_use = test;
        to_use[1] = '-'; to_use[3] = '-'; to_use[4] = '-';
        REQUIRE(to_use != test);
        REQUIRE(sp.to_string(enc.next_minimizer()) == to_use);
    }
    SECTION("small_genome") {
        v = {1, 2, 18};
        while(v.size() < 30) v.push_back(0);
        Spacer sp(31, 31, v);
        EncType enc(nullptr, 34, sp, nullptr, true);
        REQUIRE(enc.canonicalize());
        gzFile fp(gzopen("test/small_genome.fa", "rb"));
        kseq_t *ks(kseq_init(fp));
        kseq_read(ks);
        enc.assign(ks);
        REQUIRE(!enc.has_next_kmer());
        kseq_read(ks);
        enc.assign(ks);
        REQUIRE(enc.has_next_kmer());
        kseq_destroy(ks);
        gzclose(fp);
    }
    SECTION("qmap") {
        for(auto window_size: {32, 55, 100, 300, 500}) {
            std::vector<uint64_t> minimizers;
            gzFile fp(gzopen("test/phix.fa", "rb"));
            kseq_t *ks(kseq_init(fp));
            kseq_read(ks);
            v = {1, 2, 18};
            size_t vsum(0); for(auto i: v) vsum += i;
            while(v.size() < 30) v.push_back(0);
            Spacer sp(31, window_size, v);
            window_size = std::max(window_size, (int)sp.c_);
            EncType enc(ks->seq.s, ks->seq.l, sp, nullptr, true);
            REQUIRE(enc.canonicalize());
            uint64_t km, n(0);
            while(enc.has_next_kmer()) {
                ++n;
                if((km = enc.next_minimizer()) != BF)
                    minimizers.push_back(km);
            }
            REQUIRE(minimizers.size() == ks->seq.l - window_size + 1);
            gzclose(fp);
            kseq_destroy(ks);
        }
    }
    SECTION("space_case", "[no ambiguous bases]") {
        Spacer sp(31, 31);
        Encoder<score::Lex> enc(sp, true);
        gzFile fp(gzopen("test/phix.fa", "rb"));
        LOG_DEBUG("Opened fp\n");
        if(fp == nullptr) RUNTIME_ERROR("ZOMG fp is null for file at test/phix.fa");
        kseq_t *ks(kseq_init(fp));
        LOG_DEBUG("Opened kseq\n");
        if(ks == nullptr) RUNTIME_ERROR("ks is null for fp");
        std::unordered_set<u64> kmers, okmers;
        u64 k(BF);
        while(kseq_read(ks) >= 0) {
            enc.assign(ks);
            while(enc.has_next_kmer())
                if((k = enc.next_kmer()) != BF) kmers.insert(k);
        }
        LOG_DEBUG("Filled kmers\n");
        kseq_rewind(ks);
        gzrewind(fp);
        while(kseq_read(ks) >= 0) {
            LOG_DEBUG("reading seq\n");
            enc.assign(ks);
            while(enc.has_next_kmer()) {
                if((k = enc.next_kmer()) != BF) {
                    okmers.insert(k);
                }
            }
        }
        size_t olap(0);
        for(const auto k: kmers) olap += (okmers.find(k) != okmers.end());
        LOG_DEBUG("Size of each: kmers %zu, okmers %zu\n", kmers.size(), okmers.size());
        REQUIRE(kmers.size() == okmers.size());
        REQUIRE(olap == kmers.size());
        REQUIRE(kmers.size() == 5356);
        gzclose(fp);
        kseq_destroy(ks);
    }
    SECTION("space_case_jump") {
        for(const char *SPACED_FN: {"test/phix.fa", "test/ec/GCF_000007445.1_ASM744v1_genomic.fna.gz"}) {
            Spacer sp(31, 31);
            Encoder<score::Entropy> enc(sp, true);
            gzFile fp(gzopen(SPACED_FN, "rb"));
            if(fp == nullptr) RUNTIME_ERROR("Could not open file at "s + SPACED_FN);
            kseq_t *ks(kseq_init(fp));
            std::unordered_set<u64> kmers, okmers;
            u64 k(BF);
            while(kseq_read(ks) >= 0) {
                enc.assign(ks);
                while(enc.has_next_kmer())
                    if((k = enc.next_kmer()) != BF)
                        kmers.insert(canonical_representation(k, 31));
            }
            kseq_rewind(ks);
            gzrewind(fp);
            enc.set_canonicalize(false);
            enc.for_each([&](auto km) {okmers.insert(canonical_representation(km, 31));}, ks);
            size_t olap(0);
            for(const auto k: kmers) olap += (okmers.find(k) != okmers.end());
            LOG_DEBUG("Size of each: kmers %zu, okmers %zu, olap %zu\n", kmers.size(), okmers.size(), olap);
#if 0
            if(olap != kmers.size()) {
                if(okmers.size() != olap) RUNTIME_ERROR("I can't even try to help you this is so broken.");
                LOG_DEBUG("Number of missing kmers: %zu\n", size_t(kmers.size() - olap));
                for(const auto kmer: kmers)
                    if(okmers.find(kmer) == okmers.end())
                        std::fprintf(stderr, "Kmer: %s. RC: %s\n", sp.to_string(kmer).data(), sp.to_string(reverse_complement(kmer, 31)).data());
                LOG_DEBUG("Number of missing kmers: %zu\n", size_t(kmers.size() - olap));
            }
#endif
            REQUIRE(olap == kmers.size());
            REQUIRE(kmers.size() == okmers.size());
            gzclose(fp);
            kseq_destroy(ks);
        }
    }
}
TEST_CASE("rollin") {
    RollingHasher<__uint128_t, CyclicHash<__uint128_t>> enc(100);
    RollingHasher<__uint128_t, KarpRabinHashBits<__uint128_t>> enc2(100, true);
    gzFile fp = gzopen("test/phix.fa", "rb");
    if(!fp) throw "a party!";
    __uint128_t total_hash = 0;
    enc.for_each([&total_hash](auto x) {total_hash ^= x;}, fp);
    std::fprintf(stderr, "total hash sum: %zu/%zu\n", size_t(total_hash>>64), size_t(total_hash));
    gzrewind(fp);
    enc2.for_each([&total_hash](auto x) {total_hash ^= x;}, fp);
    std::fprintf(stderr, "total hash sum for canon: %zu/%zu\n", size_t(total_hash>>64), size_t(total_hash));
    gzrewind(fp);
    size_t xor_red = 0;
    Encoder<> enc3(31);
    enc3.for_each_hash([&](uint64_t x) {xor_red |= x;});
    gzclose(fp);
}
TEST_CASE("entmin") {
    Spacer sp(31, 60);
    Encoder<score::Entropy> enc(Spacer(31, 60));
    std::unordered_set<u64> kmers;
    enc.for_each([&](const u64 val) {kmers.insert(val);}, "test/phix.fa");
    REQUIRE(kmers.size() < 5353);
    std::unordered_set<u64> kmers2;
    LOG_INFO("kmers size: %zu\n", kmers.size());
    {
        Encoder<score::Entropy> enc2(Spacer(31, 31));
        enc2.for_each([&](const u64 &val) {kmers2.insert(val);}, "test/phix.fa");
    }
    LOG_INFO("kmers2 size: %zu\n", kmers2.size());
}
