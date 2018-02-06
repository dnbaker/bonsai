#include "test/catch.hpp"
#include "lib/spacer.h"
#include "lib/encoder.h"
#include <algorithm>
#include <numeric>

using namespace emp;
using EncType = Encoder<score::Lex>;

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
        EncType enc(&test[0], 34, sp, nullptr);
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
        EncType enc(&test[0], 34, sp, nullptr);
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
        EncType enc(nullptr, 34, sp, nullptr);
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
            EncType enc(ks->seq.s, ks->seq.l, sp, nullptr);
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
        Encoder<score::Entropy> enc(sp);
        gzFile fp(gzopen("test/phix.fa", "rb"));
        if(fp == nullptr) throw std::runtime_error("ZOMG fp is null for file at test/phix.fa");
        kseq_t *ks(kseq_init(fp));
        if(ks == nullptr) throw std::runtime_error("ks is null for fp");
        std::unordered_set<u64> kmers, okmers;
        u64 k(BF);
        while(kseq_read(ks) >= 0) {
            enc.assign(ks);
            while(enc.has_next_kmer()) {
                if((k = enc.next_unspaced_kmer(k)) != BF)
                    kmers.insert(k);
            }
        }
        kseq_rewind(ks);
        gzrewind(fp);
        while(kseq_read(ks) >= 0) {
            enc.assign(ks);
            while(enc.has_next_kmer()) {
                if((k = enc.next_minimizer()) != BF)
                    okmers.insert(k);
            }
        }
        LOG_DEBUG("Size of each: kmers %zu, okmers %zu\n", kmers.size(), okmers.size());
        size_t olap(0);
        for(const auto k: kmers) olap += (okmers.find(k) != okmers.end());
        REQUIRE(kmers.size() == okmers.size());
        REQUIRE(olap == kmers.size());
        gzclose(fp);
        kseq_destroy(ks);
    }
    SECTION("space_case_jump") {
        const char *SPACED_FN = "ec/GCF_000007445.1_ASM744v1_genomic.fna.gz";
        Spacer sp(31, 31);
        Encoder<score::Entropy> enc(sp);
        gzFile fp(gzopen(SPACED_FN, "rb"));
        kseq_t *ks(kseq_init(fp));
        std::unordered_set<u64> kmers, okmers;
        u64 k(BF);
        while(kseq_read(ks) >= 0) {
            enc.assign(ks);
            while(enc.has_next_kmer()) {
                if((k = enc.next_unspaced_kmer(k)) != BF)
                    kmers.insert(k);
            }
        }
        kseq_rewind(ks);
        gzrewind(fp);
        while(kseq_read(ks) >= 0) {
            LOG_DEBUG("Reading a seq from the file.\n");
            enc.assign(ks);
            while(enc.has_next_kmer()) {
                if((k = enc.next_minimizer()) != BF)
                    okmers.insert(k);
            }
        }
        LOG_DEBUG("Size of each: kmers %zu, okmers %zu\n", kmers.size(), okmers.size());
        size_t olap(0);
        for(const auto k: kmers) olap += (okmers.find(k) != okmers.end());
        REQUIRE(kmers.size() == okmers.size());
        REQUIRE(olap == kmers.size());
        gzclose(fp);
        kseq_destroy(ks);
    }
}
