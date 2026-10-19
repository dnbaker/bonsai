// Key tiers through the REAL build path (lca_map over genome files, a seq2tax
// file and a taxonomy), then classify with the khash and the compact database:
//  - k = 31: exact 64-bit packed k-mers (KmerKeys::EXACT64), 62-bit keys;
//  - k = 40: 61-bit rolling-hash keys (KmerKeys::HASHED64; formerly rejected).
// For each k: the database keys span the expected width, kraken-style output
// is byte-identical between khash and compact databases, and error-free reads
// classify to their own genome's taxon (which fails if build and classify
// ever disagree on keys).
#include "test/catch.hpp"

#include "classifier.h"
#include "feature_min.h"
#include "kmermap.h"

#include <cstdio>
#include <map>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using namespace bns;

namespace {

std::string random_dna(size_t n, std::mt19937_64 &rng) {
    static const char b[] = "ACGT";
    std::string s(n, 'A');
    for(auto &c: s) c = b[rng() & 3];
    return s;
}
std::string slurp(std::FILE *fp) {
    std::fflush(fp);
    std::rewind(fp);
    std::string ret;
    char buf[1 << 14];
    size_t n;
    while((n = std::fread(buf, 1, sizeof(buf), fp)) > 0) ret.append(buf, n);
    return ret;
}
void put(const std::string &path, const std::string &content) {
    std::FILE *fp = std::fopen(path.data(), "w");
    REQUIRE(fp != nullptr);
    std::fputs(content.data(), fp);
    std::fclose(fp);
}

} // namespace

TEST_CASE("KmerKeys tiers by alphabet capacity") {
    REQUIRE(KmerKeys<>::tier_for(32) == KmerKeys<>::EXACT64);
    REQUIRE(KmerKeys<>::tier_for(33) == KmerKeys<>::HASHED64);
    REQUIRE(KmerKeys<>::key_bits_for(31) == 62);
    REQUIRE(KmerKeys<>::key_bits_for(32) == 64);
    REQUIRE(KmerKeys<>::key_bits_for(40) == 61);
    REQUIRE(KmerKeys<>::key_bits_for(40, DNA, 40) == 40);
    // bonsai's PROTEIN packs 8 bits/residue (8 per word); PROTEIN20 is the
    // base-20 integer encoding (14 per word: 20^14 < 2^61)
    REQUIRE(KmerKeys<>::tier_for(8, PROTEIN) == KmerKeys<>::EXACT64);
    REQUIRE(KmerKeys<>::tier_for(9, PROTEIN) == KmerKeys<>::HASHED64);
    REQUIRE(KmerKeys<>::tier_for(14, PROTEIN20) == KmerKeys<>::EXACT64);
    REQUIRE(KmerKeys<>::tier_for(15, PROTEIN20) == KmerKeys<>::HASHED64);
    REQUIRE(KmerKeys<>::key_bits_for(14, PROTEIN20) == 61);   // bit width of 20^14 - 1
    using K128 = KmerKeys<score::Lex, u128>;
    REQUIRE(K128::tier_for(40) == K128::EXACT128);
    REQUIRE(K128::tier_for(64) == K128::EXACT128);
    REQUIRE(K128::tier_for(65) == K128::HASHED128);
    REQUIRE(K128::key_bits_for(40) == 80);
    REQUIRE(K128::key_bits_for(64) == 128);
    REQUIRE(K128::key_bits_for(80) == 122);
}

TEST_CASE("lca_map + classify for exact (k = 31) and hashed (k = 40) keys") {
    char tmpl[] = "/tmp/bns_hk_XXXXXX";
    const std::string dir(mkdtemp(tmpl));
    std::mt19937_64 rng(13);
    std::vector<std::string> genomes, fns;
    std::string seq2tax;
    for(int g = 0; g < 4; ++g) {
        genomes.push_back(random_dna(40000, rng));
        const std::string fn = dir + "/g" + std::to_string(g) + ".fa";
        put(fn, ">g" + std::to_string(g) + " synthetic\n" + genomes.back() + "\n");
        fns.push_back(fn);
        seq2tax += "g" + std::to_string(g) + "\t" + std::to_string(2 + g) + "\n";
    }
    put(dir + "/seq2tax.txt", seq2tax);
    put(dir + "/tax.txt", "2\t|\t1\n3\t|\t1\n4\t|\t1\n5\t|\t1\n");
    khash_t(p) *taxmap(build_parent_map((dir + "/tax.txt").data()));
    // reads: error-free substrings (for the taxon check), plus N-containing and random reads
    const std::string fq = dir + "/r.fq";
    std::map<std::string, tax_t> truth;
    {
        std::FILE *a = std::fopen(fq.data(), "w");
        for(int r = 0; r < 1200; ++r) {
            std::string s;
            const std::string name = "r" + std::to_string(r);
            if(r % 10 == 9) s = random_dna(150, rng);
            else {
                const int g = int(rng() % 4);
                s = genomes[g].substr(rng() % (genomes[g].size() - 150), 150);
                if(r % 10 == 8) s[75] = 'N';
                else truth[name] = tax_t(2 + g);
            }
            std::fprintf(a, "@%s\n%s\n+\n%s\n", name.data(), s.data(), std::string(s.size(), 'I').data());
        }
        std::fclose(a);
    }
    for(const unsigned k: {31u, 40u}) {
        Spacer sp(k, k);
        khash_t(c) *db(lca_map<score::Lex>(fns, taxmap, (dir + "/seq2tax.txt").data(), sp, 2, true, 1 << 16));
        REQUIRE(kh_size(db) > 100000);
        const unsigned kb = KmerKeys<>::key_bits_for(k);
        REQUIRE(CompactKmerMap::key_bits_needed(db) <= kb);
        CompactKmerMap cm(db, CompactKmerMap::DEFAULT_MAX_LOAD, kb);
        REQUIRE(cm.key_bits() == kb);
        ClassifierGeneric<score::Lex> ck(db, spvec_t{}, k, k, 2, true, false, true, true);
        ClassifierGeneric<score::Lex, CompactKmerMap> cc(&cm, spvec_t{}, k, k, 2, true, false, true, true);
        REQUIRE(ck.enc_.tier() == (k <= 32 ? KmerKeys<>::EXACT64: KmerKeys<>::HASHED64));
        std::FILE *ok = std::tmpfile(), *oc = std::tmpfile();
        process_dataset(ck, taxmap, fq.data(), nullptr, ok, 256, 32);
        process_dataset(cc, taxmap, fq.data(), nullptr, oc, 256, 32);
        const std::string a = slurp(ok), b = slurp(oc);
        std::fclose(ok); std::fclose(oc);
        REQUIRE(a == b);
        size_t right = 0, total = 0;
        std::istringstream is(a);
        std::string line;
        while(std::getline(is, line)) {
            std::istringstream ls(line);
            std::string cu, name;
            tax_t taxon = 0;
            ls >> cu >> name >> taxon;
            const auto t = truth.find(name);
            if(t == truth.end()) continue;
            ++total;
            right += taxon == t->second;
        }
        REQUIRE(total == truth.size());
        REQUIRE(right == total);           // random genomes: every error-free read is unique to its genome
        // direct compact build (no khash), from a good estimate and from a tiny
        // one that forces repeated growth: identical content to the khash build
        for(const size_t est: {size_t(kh_size(db)), size_t(1000)}) {
            auto dm = compact_lca_map<score::Lex>(fns, taxmap, (dir + "/seq2tax.txt").data(), sp, 2, true, est, kb);
            REQUIRE(dm->size() == kh_size(db));
            REQUIRE(dm->key_bits() == kb);
            REQUIRE(dm->value_bits() == CompactKmerMap::bits_for(dm->num_taxa()));
            size_t bad = 0;
            for(khint_t ki = kh_begin(db); ki != kh_end(db); ++ki)
                if(kh_exist(db, ki)) bad += dm->get(kh_key(db, ki)) != kh_val(db, ki);
            REQUIRE(bad == 0);
            ClassifierGeneric<score::Lex, CompactKmerMap> cd(dm.get(), spvec_t{}, k, k, 2, true, false, true, true);
            std::FILE *od = std::tmpfile();
            process_dataset(cd, taxmap, fq.data(), nullptr, od, 256, 32);
            REQUIRE(slurp(od) == a);
            std::fclose(od);
            std::fprintf(stderr, "  direct compact build from estimate %zu: %zu k-mers, %zu bytes, identical to khash build\n",
                         est, dm->size(), dm->bytes_used());
        }
        std::fprintf(stderr, "k=%u (%s keys, %u bits): %zu k-mers; %zu/%zu error-free reads to their own taxon; khash == compact\n",
                     k, k <= 32 ? "exact": "hashed", kb, size_t(kh_size(db)), right, total);
        if(k == 40) {
            // EXACT 128-bit keys (32 < k <= 64): compact build only. Its classify
            // output must equal the 61-bit-hashed databases' (identical k-mer
            // sets -> identical LCAs, absent hash collisions at this size), and
            // a forced-growth build must equal a well-sized one, key for key.
            const unsigned kb128 = KmerKeys<score::Lex, u128>::key_bits_for(k);
            auto m128 = compact_lca_map<score::Lex, u128>(fns, taxmap, (dir + "/seq2tax.txt").data(), sp, 2, true, kh_size(db), kb128);
            auto g128 = compact_lca_map<score::Lex, u128>(fns, taxmap, (dir + "/seq2tax.txt").data(), sp, 2, true, 1000, kb128);
            REQUIRE(m128->key_bits() == 80);
            REQUIRE(m128->size() == kh_size(db));          // as many exact k-mers as hashed keys: no hash merges here
            REQUIRE(g128->size() == m128->size());
            size_t diff = 0;
            m128->for_each([&](u128 km, tax_t t) {diff += g128->get(km) != t;});
            REQUIRE(diff == 0);
            const std::string cdb = dir + "/db128.cdb";
            CompactDatabase::write(cdb.data(), k, k, spvec_t(k - 1, 0), *m128);
            CompactDatabase cd(cdb.data());
            REQUIRE(cd.wide());
            REQUIRE(cd.db128_->size() == m128->size());
            ClassifierGeneric<score::Lex, CompactKmerMap128> c128(cd.db128_.get(), spvec_t{}, k, k, 2, true, false, true, true);
            REQUIRE(c128.enc_.tier() == KmerKeys<score::Lex, u128>::EXACT128);
            std::FILE *o128 = std::tmpfile();
            process_dataset(c128, taxmap, fq.data(), nullptr, o128, 256, 32);
            REQUIRE(slurp(o128) == a);
            std::fclose(o128);
            std::remove(cdb.data());
            std::fprintf(stderr, "  k=40 exact 128-bit keys (80 bits): %zu k-mers, %zu bytes; classify == hashed-key databases\n",
                         m128->size(), m128->bytes_used());
        }
        kh_destroy(c, db);
    }
    // SPACED seeds and minimizer WINDOWS beyond k = 32 (32 < k <= 64): the
    // 61-bit keys fold the exact 128-bit spaced k-mer, so khash, compact
    // (61-bit) and exact 128-bit databases select the same keys and classify
    // identically; error-free reads still reach their own taxon.
    {
        const unsigned k = 40;
        spvec_t gaps(k - 1, 0);
        for(unsigned i = 4; i < k - 1; i += 6) gaps[i] = 1;      // a gap after every 6th position: span > k
        for(const unsigned w: {k + 6u, k + 16u}) {
            Spacer sp(k, w, gaps);
            REQUIRE(!sp.unspaced());
            khash_t(c) *db(lca_map<score::Lex>(fns, taxmap, (dir + "/seq2tax.txt").data(), sp, 2, true, 1 << 16));
            auto m128 = compact_lca_map<score::Lex, u128>(fns, taxmap, (dir + "/seq2tax.txt").data(), sp, 2, true, kh_size(db),
                                                          KmerKeys<score::Lex, u128>::key_bits_for(k));
            REQUIRE(m128->size() == kh_size(db));
            CompactKmerMap cm(db, CompactKmerMap::DEFAULT_MAX_LOAD, KmerKeys<>::key_bits_for(k));
            spvec_t sv(gaps);
            ClassifierGeneric<score::Lex> ck(db, sv, k, w, 2, true, false, true, true);
            ClassifierGeneric<score::Lex, CompactKmerMap> cc(&cm, sv, k, w, 2, true, false, true, true);
            ClassifierGeneric<score::Lex, CompactKmerMap128> c8(m128.get(), sv, k, w, 2, true, false, true, true);
            REQUIRE(ck.enc_.tier() == KmerKeys<>::HASHED64);
            REQUIRE(c8.enc_.tier() == KmerKeys<score::Lex, u128>::EXACT128);
            std::string out[3];
            int ix = 0;
            auto run = [&](const auto &cl) {
                std::FILE *o = std::tmpfile();
                process_dataset(cl, taxmap, fq.data(), nullptr, o, 256, 32);
                out[ix++] = slurp(o);
                std::fclose(o);
            };
            run(ck); run(cc); run(c8);
            REQUIRE(out[0] == out[1]);
            REQUIRE(out[0] == out[2]);
            size_t right = 0, total = 0;
            std::istringstream is(out[0]);
            std::string line;
            while(std::getline(is, line)) {
                std::istringstream ls(line);
                std::string cu, name;
                tax_t taxon = 0;
                ls >> cu >> name >> taxon;
                const auto t = truth.find(name);
                if(t == truth.end()) continue;
                ++total; right += taxon == t->second;
            }
            REQUIRE(total == truth.size());
            REQUIRE(right == total);
            std::fprintf(stderr, "k=40 spaced (span %u), window %u: %zu keys; 61-bit khash == 61-bit compact == exact 128-bit; %zu/%zu reads to own taxon\n",
                         sp.c(), w, size_t(kh_size(db)), right, total);
            kh_destroy(c, db);
        }
    }
    kh_destroy(p, taxmap);
}
