#include "test/catch.hpp"
#include "spacer.h"
#include "encoder.h"
#include <algorithm>
#include <numeric>
#include <glob.h>

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
        REQUIRE(!enc.canonicalize()); // Spaced can't be canon
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
        REQUIRE(!enc.canonicalize()); // Spaced kmers cannot be canonicalized (unless they're palindromic, but I haven't enabled that edge case.)
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
        REQUIRE(!enc.canonicalize());
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
            REQUIRE(!enc.canonicalize());
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
        std::unordered_map<u64, u32> kmers, okmers;
        u64 k(BF);
        while(kseq_read(ks) >= 0) {
            enc.assign(ks);
            while(enc.has_next_kmer())
                if((k = enc.next_kmer()) != BF) ++kmers[k];
        }
        kseq_destroy(ks);
        gzclose(fp);
        fp = gzopen("test/phix.fa", "rb");
        ks = kseq_init(fp);
        LOG_DEBUG("Filled kmers\n");
        while(kseq_read(ks) >= 0) {
            enc.assign(ks);
            while(enc.has_next_kmer()) {
                if((k = enc.next_kmer()) != BF)
                    ++okmers[k];
            }
        }
        size_t olap(0);
        for(const auto &k: kmers) olap += (okmers.find(k.first) != okmers.end());
        LOG_DEBUG("Size of each: kmers %zu, okmers %zu\n", kmers.size(), okmers.size());
        REQUIRE(kmers.size() == okmers.size());
        REQUIRE(olap == kmers.size());
        REQUIRE(kmers.size() == 5356);
        gzclose(fp);
        kseq_destroy(ks);
    }
}
TEST_CASE("space_case_jump") {
    for(const char *SPACED_FN: {"test/phix.fa", "test/ec/GCF_000007445.1_ASM744v1_genomic.fna.gz"}) {
        Spacer sp(31, 31);
        Encoder<score::Entropy> enc(sp);
        gzFile fp(gzopen(SPACED_FN, "rb"));
        if(fp == nullptr) RUNTIME_ERROR("Could not open file at "s + SPACED_FN);
        kseq_t *ks(kseq_init(fp));
        std::unordered_set<u64> kmers, okmers;
        //u64 k(BF);
        enc.for_each_canon([&](auto km) {kmers.insert(km);}, ks);
        kseq_destroy(ks);
        gzclose(fp);
        fp = gzopen(SPACED_FN, "rb");
        ks = kseq_init(fp);
        //enc.canonicalize(false);
        enc.for_each_uncanon([&](auto km) {okmers.insert(canonical_representation(km, 31));}, ks);
        size_t olap(0);
        for(const auto k: kmers) olap += (okmers.find(k) != okmers.end());
        LOG_DEBUG("Size of each: kmers %zu, okmers %zu, olap %zu\n", kmers.size(), okmers.size(), olap);
        REQUIRE(olap == kmers.size());
        REQUIRE(kmers.size() == okmers.size());
        gzclose(fp);
        kseq_destroy(ks);
    }
}
TEST_CASE("rollin_spaced") {
    RollingHasher<__uint128_t, CyclicHash<__uint128_t>> enc(100, /*canon = */false, /*alphabet=*/bns::PROTEIN, /*wsz=*/200);
    __uint128_t total_hash = 0;
    enc.for_each_hash([&total_hash](auto) {++total_hash;}, "test/phix.fa");
    REQUIRE(uint64_t(total_hash) == uint64_t(5386 - 200 + 1));
    size_t xor_red = 0;
    Encoder<> enc3(31);
    enc3.for_each_hash([&](uint64_t x) {xor_red |= x;}, "test/phix.fa");
    Encoder<> enc4(147);
    enc3.for_each_hash([&](uint64_t x) {xor_red |= x;}, "test/phix.fa");
}
TEST_CASE("rollin") {
    RollingHasher<__uint128_t, CyclicHash<__uint128_t>> enc(100);
    //RollingHasher<__uint128_t, KarpRabinHashBits<__uint128_t>> enc2(100, true);
    gzFile fp = gzopen("test/phix.fa", "rb");
    if(!fp) throw "a party!";
    __uint128_t total_hash = 0;
    enc.for_each_hash([&total_hash](auto x) {total_hash ^= x;}, fp);
    std::fprintf(stderr, "total hash sum: %zu/%zu\n", size_t(total_hash>>64), size_t(total_hash));
    gzrewind(fp);
    //enc2.for_each_canon([&total_hash](auto x) {total_hash ^= x;}, fp);
    std::fprintf(stderr, "total hash sum for canon: %zu/%zu\n", size_t(total_hash>>64), size_t(total_hash));
    gzrewind(fp);
    size_t xor_red = 0;
    Encoder<> enc3(31);
    enc3.for_each_hash([&](uint64_t x) {xor_red |= x;});
    Encoder<> enc4(147);
    enc3.for_each_hash([&](uint64_t x) {xor_red |= x;});
    gzclose(fp);
}
TEST_CASE("rhs") {
    RollingHasherSet<uint64_t> rhs(std::vector<int>{16,32,64,128});
    std::FILE *ofp = std::fopen("/dev/null", "wb");
    rhs.for_each_hash([ofp](uint64_t kmer, size_t hashnum) {std::fprintf(ofp, "%" PRIu64 ": %zu\n", kmer, hashnum);},
                      "test/phix.fa");
    std::fclose(ofp);
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
TEST_CASE("parseasprot") {
    Spacer sp(31, 31);
    for(const InputType rv: {bns::DNA, bns::PROTEIN, bns::PROTEIN20, bns::PROTEIN_3BIT, bns::PROTEIN_14, bns::PROTEIN_6, bns::DNA2}) {
        Encoder<score::Lex> enc(sp, true);
        enc.hashtype(rv);
        gzFile fp(gzopen("test/phix.fa", "rb"));
        LOG_DEBUG("Opened fp\n");
        if(fp == nullptr) RUNTIME_ERROR("ZOMG fp is null for file at test/phix.fa");
        kseq_t *ks(kseq_init(fp));
        LOG_INFO("Opened kseq\n");
        if(ks == nullptr) RUNTIME_ERROR("ks is null for fp");
        std::unordered_map<u64, u32> kmers, okmers;
        u64 k(BF);
        while(kseq_read(ks) >= 0) {
            enc.assign(ks);
            while(enc.has_next_kmer())
                if((k = enc.next_kmer()) != BF) ++kmers[k];
        }
        LOG_INFO("Filled kmers\n");
        kseq_destroy(ks); gzclose(fp);
        fp = gzopen("test/phix.fa", "rb");
        ks = kseq_init(fp);
        while(kseq_read(ks) >= 0) {
            enc.assign(ks);
            while(enc.has_next_kmer()) {
                if((k = enc.next_kmer()) != BF) { ++okmers[k];}
            }
        }
        size_t olap(0);
        for(const auto &pair: kmers) olap += (okmers.find(pair.first) != okmers.end());
        LOG_INFO("Size of each: kmers %zu, okmers %zu\n", kmers.size(), okmers.size());
        gzclose(fp);
        kseq_destroy(ks);
        std::fprintf(stderr, "Finished for %s\n", bns::to_string(rv).data());
    }
}
std::string gz2xz(std::string x) {
    x[x.size() - 2] = 'x';
    return x;
}
std::string gz2bz(std::string x) {
    x = x.substr(0, x.find_last_of('.'));
    x = x + ".bz2";
    return x;
}
std::string gz2zst(std::string x) {
    x = x.substr(0, x.find_last_of('.'));
    x = x + ".zst";
    return x;
}

TEST_CASE("xzparse") {
    Spacer sp(31, 71);
    glob_t glo;
    std::memset(&glo, 0, sizeof(glo));
    if(::glob("test/ec/*fna.gz", GLOB_TILDE, NULL, &glo)) {
        throw std::runtime_error("Glob failed");
    }
    std::vector<std::string> xzs, gzs, bzs, zzs;
    for(size_t i = 0; i < glo.gl_pathc; ++i) {
        gzs.emplace_back(glo.gl_pathv[i]);
        xzs.emplace_back(gz2xz(gzs.back()));
        bzs.emplace_back(gz2bz(gzs.back()));
        zzs.emplace_back(gz2zst(gzs.back()));
    }
    std::FILE *tmp = ::popen("zstd -h 2>/dev/null", "r");
    const bool haszstd = tmp && !::pclose(tmp);
    OMP_PFOR
    for(size_t i = 0; i < glo.gl_pathc; ++i) {
        std::FILE *fp;
        std::string cmd;
        if(!bns::isfile(xzs[i])) {
            cmd = std::string("ls ") + xzs[i] + " &>/dev/null || " + "gzip -dc " + gzs[i] + " | xz > " + xzs[i];
            if(!(fp = ::popen(cmd.data(), "r"))) throw 1;
            ::pclose(fp);
        }
        if(!bns::isfile(bzs[i])) {
            cmd = std::string("ls ") + bzs[i] + " &>/dev/null || " + "gzip -dc " + gzs[i] + " | bzip2 > " + bzs[i];
            if((fp = ::popen(cmd.data(), "r")) == nullptr) throw 1;
            ::pclose(fp);
        }
        if(haszstd && !bns::isfile(zzs[i])) {
            cmd = std::string("ls ") + zzs[i] + " &>/dev/null || " + "gzip -dc " + gzs[i] + " | zstd > " + zzs[i];
            if((fp = ::popen(cmd.data(), "r")) == nullptr) throw 1;
            ::pclose(fp);
        }
    }
    OMP_PFOR
    for(size_t i = 0; i < glo.gl_pathc; ++i) {
        uint64_t lhv = 0, rhv = 0, bhv = 0, zhv = 0;
        {
            Encoder<> enc(sp);
            enc.for_each([&lhv](auto x) {lhv ^= x;}, gzs[i].data());
        }
        {
            Encoder<> enc(sp);
            enc.for_each([&rhv](auto x) {rhv ^= x;}, xzs[i].data());
        }
        {
            Encoder<> enc(sp);
            enc.for_each([&bhv](auto x) {bhv ^= x;}, bzs[i].data());
        }
        if(haszstd) {
            Encoder<> enc(sp);
            enc.for_each([&zhv](auto x) {zhv ^= x;}, zzs[i].data());
        } else {
            zhv = lhv;
        }
        REQUIRE(lhv == rhv);
        REQUIRE(lhv == bhv);
        REQUIRE(lhv == zhv);
        lhv = rhv = bhv = zhv = 0;
        {
            RollingHasher<uint64_t> rh(31, false);
            rh.for_each([&lhv](auto x) {lhv ^= x;}, gzs[i].data());
        }
        {
            RollingHasher<uint64_t> rh(31, false);
            rh.for_each([&rhv](auto x) {rhv ^= x;}, xzs[i].data());
        }
        {
            RollingHasher<uint64_t> rh(31, false);
            rh.for_each([&bhv](auto x) {bhv ^= x;}, bzs[i].data());
        }
        if(haszstd) {
            RollingHasher<uint64_t> rh(31, false);
            rh.for_each([&zhv](auto x) {zhv ^= x;}, zzs[i].data());
        } else {
            zhv = lhv;
        }
        REQUIRE(lhv == rhv);
        REQUIRE(lhv == bhv);
        REQUIRE(lhv == zhv);
    }
    globfree(&glo);
}
