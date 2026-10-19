// CompactKmerMap (kmermap.h) must be a drop-in for the khash_t(c) k-mer -> taxon
// database:
//  - identical get/get_many results, NOT_FOUND for absent k-mers, tax 0 kept,
//    overwrite, and tax_t(-1) rejected;
//  - a CompactDatabase write/read round trip;
//  - END TO END: classification output (kraken-style, emit-all) must be
//    byte-identical with the khash and the compact database, single- and
//    paired-end, on reads drawn from the database genomes with substitutions
//    and Ns, plus random and reverse-complement reads.
#include "test/catch.hpp"

#include "classifier.h"
#include "kmermap.h"

#include <cstdio>
#include <random>
#include <unordered_set>
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
std::string revcomp(const std::string &s) {
    std::string r(s.rbegin(), s.rend());
    for(auto &c: r) c = c == 'A' ? 'T': c == 'C' ? 'G': c == 'G' ? 'C': c == 'T' ? 'A': c;
    return r;
}
std::string tmpdir() {
    char tmpl[] = "/tmp/bns_ckm_XXXXXX";
    const char *ret(mkdtemp(tmpl));
    REQUIRE(ret != nullptr);
    return ret;
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

} // namespace

TEST_CASE("CompactKmerMap matches khash_t(c) lookups") {
    std::mt19937_64 rng(5);
    khash_t(c) *km(kh_init(c));
    int khr;
    for(int i = 0; i < 300000; ++i) {
        const u64 k = rng() % 250000;              // repeats -> overwrites
        const tax_t t = (i % 1000 == 0) ? 0: tax_t(rng() % 3300000u);
        const khint_t ki(kh_put(c, km, k, &khr));
        kh_val(km, ki) = t;
    }
    CompactKmerMap cm(km);
    REQUIRE(cm.size() == kh_size(km));
    std::vector<u64> qs;
    for(int i = 0; i < 200000; ++i) qs.push_back((i & 1) ? rng(): u64(rng() % 250000));
    std::vector<tax_t> got(qs.size());
    std::vector<uint64_t> scratch;
    cm.get_many(qs.begin(), qs.end(), got.begin(), scratch);
    for(size_t i = 0; i < qs.size(); ++i) {
        const khint_t ki = kh_get(c, km, qs[i]);
        const tax_t want = ki == kh_end(km) ? CompactKmerMap::NOT_FOUND: kh_val(km, ki);
        REQUIRE(cm.get(qs[i]) == want);
        REQUIRE(got[i] == want);
    }
    size_t n = 0, zeros = 0;
    cm.for_each([&](u64 k, tax_t t) {
        ++n; zeros += t == 0;
        REQUIRE(kh_val(km, kh_get(c, km, k)) == t);
    });
    REQUIRE(n == kh_size(km));
    REQUIRE(zeros > 0);
    REQUIRE_THROWS_AS(cm.set(1, tax_t(-1)), std::invalid_argument);
    // dense ids: the minimum width for the distinct taxa
    {
        std::unordered_set<tax_t> distinct;
        for(khint_t ki = kh_begin(km); ki != kh_end(km); ++ki) if(kh_exist(km, ki)) distinct.insert(kh_val(km, ki));
        REQUIRE(cm.num_taxa() == distinct.size());
        REQUIRE(cm.value_bits() == CompactKmerMap::bits_for(distinct.size()));
        REQUIRE(((uint64_t(1) << (cm.value_bits() - 1)) - 1) < distinct.size());
    }
    // overwrite through set() after finalize() (index rebuilt on demand)
    {
        const u64 k0 = qs[0];
        cm.finalize();
        tax_t some = 0;
        cm.for_each([&](u64, tax_t t) {some = t;});
        cm.set(k0, some);
        REQUIRE(cm.get(k0) == some);
    }
    std::fprintf(stderr, "CompactKmerMap: %zu k-mers, %zu bytes vs khash_t(c) %zu bytes (%.1f%%)\n",
                 cm.size(), cm.bytes_used(), khash_bytes(km), 100. * cm.bytes_used() / khash_bytes(km));
    kh_destroy(c, km);
}

TEST_CASE("CompactKmerMap from khash with large custom taxids (hash-map index path)") {
    khash_t(c) *km(kh_init(c));
    int khr;
    std::mt19937_64 rng(77);
    for(int i = 0; i < 50000; ++i) {
        const khint_t ki(kh_put(c, km, rng(), &khr));
        kh_val(km, ki) = tax_t(2000000000u + rng() % 4000u);   // above the direct-array cap
    }
    CompactKmerMap cm(km);
    REQUIRE(cm.num_taxa() <= 4000);
    REQUIRE(cm.value_bits() == CompactKmerMap::bits_for(cm.num_taxa()));
    for(khint_t ki = kh_begin(km); ki != kh_end(km); ++ki)
        if(kh_exist(km, ki)) REQUIRE(cm.get(kh_key(km, ki)) == kh_val(km, ki));
    kh_destroy(c, km);
}

TEST_CASE("CompactKmerMap incremental build: taxon capacity and widths") {
    CompactKmerMap m(1000, 5);                       // room for 5 distinct taxa -> 3-bit values
    REQUIRE(m.value_bits() == 3);
    for(u64 k = 0; k < 500; ++k) m.set(k * 977, tax_t(100 + k % 5));
    REQUIRE(m.num_taxa() == 5);
    REQUIRE_THROWS_AS(m.set(123456789, tax_t(999)), std::length_error);
    for(u64 k = 0; k < 500; ++k) REQUIRE(m.get(k * 977) == tax_t(100 + k % 5));
    REQUIRE(m.get(1) == CompactKmerMap::NOT_FOUND);
    // taxon 0 would be a sixth distinct taxon: rejected, old value untouched
    REQUIRE_THROWS_AS(m.set(0, 0), std::length_error);
    REQUIRE(m.get(0) == tax_t(100));
    CompactKmerMap z(100, 2);                        // tax 0 is an ordinary taxon
    z.set(5, 0); z.set(6, 9);
    REQUIRE(z.get(5) == 0);
    REQUIRE(z.get(6) == 9);
}

TEST_CASE("classify output is identical with khash and compact databases") {
    std::mt19937_64 rng(9);
    const unsigned k = 31;
    // Four "genomes" (taxa 2..5 under root 1).
    std::vector<std::string> genomes;
    for(int g = 0; g < 4; ++g) genomes.push_back(random_dna(60000, rng));
    khash_t(p) *taxmap(kh_init(p));
    int khr;
    {
        khint_t ki = kh_put(p, taxmap, 1, &khr); kh_val(taxmap, ki) = 0;
        for(tax_t t = 2; t <= 5; ++t) {ki = kh_put(p, taxmap, t, &khr); kh_val(taxmap, ki) = 1;}
    }
    // Database: every canonical k-mer of each genome, plus a shared block
    // (genome 0's first 5 kb also inserted into genome 1, so its k-mers' LCA
    // is the root 1).
    khash_t(c) *db(kh_init(c));
    ClassifierGeneric<score::Lex> builder(db, spvec_t{}, k, k, 1, true, false, true, true);
    Encoder<score::Lex> enc(builder.enc_);
    for(int g = 0; g < 4; ++g) {
        enc.for_each([&](u64 km) {
            const khint_t ki = kh_put(c, db, km, &khr);
            kh_val(db, ki) = khr == 0 && kh_val(db, ki) != tax_t(2 + g) ? 1: tax_t(2 + g);
        }, genomes[g].data(), genomes[g].size());
    }
    const std::string shared = genomes[0].substr(0, 5000);
    enc.for_each([&](u64 km) {
        const khint_t ki = kh_put(c, db, km, &khr);
        if(khr == 0 && kh_val(db, ki) != 3) kh_val(db, ki) = 1;
    }, shared.data(), shared.size());
    // Reads.
    const std::string dir = tmpdir();
    const std::string fq1 = dir + "/r1.fq", fq2 = dir + "/r2.fq";
    {
        std::FILE *a = std::fopen(fq1.data(), "w"), *b = std::fopen(fq2.data(), "w");
        REQUIRE(a); REQUIRE(b);
        for(int r = 0; r < 3000; ++r) {
            std::string s1, s2;
            if(r % 6 == 5) {s1 = random_dna(150, rng); s2 = random_dna(150, rng);}
            else {
                const auto &g = genomes[rng() % 4];
                const size_t p = rng() % (g.size() - 450);
                s1 = g.substr(p, 150); s2 = revcomp(g.substr(p + 250, 150));
                for(auto *s: {&s1, &s2})
                    for(auto &ch: *s) {
                        const auto u = rng() % 1000;
                        if(u < 20) ch = "ACGT"[rng() & 3];
                        else if(u < 23) ch = 'N';
                    }
                if(r % 6 == 4) s1 = revcomp(s1);
            }
            std::fprintf(a, "@r%d/1\n%s\n+\n%s\n", r, s1.data(), std::string(s1.size(), 'I').data());
            std::fprintf(b, "@r%d/2\n%s\n+\n%s\n", r, s2.data(), std::string(s2.size(), 'I').data());
        }
        std::fclose(a); std::fclose(b);
    }
    // Compact database through a file round trip.
    const std::string cdb = dir + "/db.cdb";
    {
        REQUIRE(CompactKmerMap::key_bits_needed(db) <= 2 * k);
        CompactKmerMap m(db, CompactKmerMap::DEFAULT_MAX_LOAD, 2 * k);   // 62-bit key domain for k = 31
        REQUIRE(m.key_bits() == 2 * k);
        // out-of-domain keys: rejected by set(), absent for get/get_many (never aliased)
        const u64 outside = (u64(1) << 62) | 12345;
        REQUIRE_THROWS_AS(m.set(outside, 2), std::invalid_argument);
        REQUIRE(m.get(outside) == CompactKmerMap::NOT_FOUND);
        std::vector<u64> qs{outside, outside & ((u64(1) << 62) - 1)};
        std::vector<tax_t> got(2);
        std::vector<uint64_t> scratch;
        m.get_many(qs.begin(), qs.end(), got.begin(), scratch);
        REQUIRE(got[0] == CompactKmerMap::NOT_FOUND);
        REQUIRE(got[1] == m.get(qs[1]));
        m.finalize();
        CompactDatabase::write(cdb.data(), k, k, spvec_t(k - 1, 0), m);
    }
    REQUIRE(CompactDatabase::is_compact(cdb.data()));
    CompactDatabase cd(cdb.data());
    REQUIRE(cd.k_ == k);
    REQUIRE(cd.db_->key_bits() == 2 * k);
    REQUIRE(cd.db_->table().cell_bits() == 2 * k - cd.db_->table().q() + cd.db_->value_bits() + 8);
    REQUIRE(cd.db_->size() == kh_size(db));
    REQUIRE(cd.db_->table().structure_check());

    ClassifierGeneric<score::Lex> ck(db, spvec_t{}, k, k, 2, true, false, true, true);
    ClassifierGeneric<score::Lex, CompactKmerMap> cc(cd.db_.get(), spvec_t{}, k, k, 2, true, false, true, true);
    for(const bool paired: {false, true}) {
        std::FILE *ok = std::tmpfile(), *oc = std::tmpfile();
        process_dataset(ck, taxmap, fq1.data(), paired ? fq2.data(): nullptr, ok, 512, 32);
        process_dataset(cc, taxmap, fq1.data(), paired ? fq2.data(): nullptr, oc, 512, 32);
        const std::string a = slurp(ok), b = slurp(oc);
        std::fclose(ok); std::fclose(oc);
        REQUIRE(a.size() > 1000);
        REQUIRE(a == b);
        // non-trivial: most reads classified, some to the root, some unclassified
        REQUIRE(a.find("\t1\t") != std::string::npos);
        std::fprintf(stderr, "classify %s: %zu bytes of output identical (khash vs compact)\n",
                     paired ? "paired": "single", a.size());
    }
    REQUIRE(ck.n_classified() == cc.n_classified());
    REQUIRE(ck.n_unclassified() == cc.n_unclassified());
    REQUIRE(ck.n_unclassified() > 0);
    kh_destroy(c, db);
    kh_destroy(p, taxmap);
    std::remove(fq1.data()); std::remove(fq2.data()); std::remove(cdb.data());
}
