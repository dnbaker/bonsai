// Kraken2-style confidence walk-up and minimum hit groups (classifier.h), with
// exact expectations, on khash and compact databases (whose batched path also
// collapses consecutive equal keys):
//  - chimera of 100 bp from genome g0 (taxon 2) + 50 bp from g1 (taxon 3) at
//    k = 31: 70 hits to 2, 20 to 3, 30 boundary misses (120 k-mers). The call
//    is 2. confidence 0.5 keeps 2 (70/120), 0.7 walks up to the root 1
//    (90/120), 0.8 leaves it unclassified;
//  - an error-free 150 bp read has 120 distinct hitting keys: classified at
//    min_hit_groups 120, unclassified at 121;
//  - defaults (0, 0) leave the output unchanged; confidence 1.0 classifies
//    every error-free read and no random read.
#include "test/catch.hpp"

#include "classifier.h"
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
template<typename CL>
std::map<std::string, tax_t> run_and_parse(const CL &cl, const khash_t(p) *taxmap, const std::string &fq, std::string *raw = nullptr) {
    std::FILE *o = std::tmpfile();
    process_dataset(cl, taxmap, fq.data(), nullptr, o, 256, 32);
    std::fflush(o); std::rewind(o);
    std::string out;
    char buf[1 << 14];
    size_t n;
    while((n = std::fread(buf, 1, sizeof(buf), o)) > 0) out.append(buf, n);
    std::fclose(o);
    if(raw) *raw = out;
    std::map<std::string, tax_t> ret;
    std::istringstream is(out);
    std::string line;
    while(std::getline(is, line)) {
        std::istringstream ls(line);
        std::string cu, name;
        tax_t taxon = 0;
        ls >> cu >> name >> taxon;
        ret[name] = taxon;
    }
    return ret;
}

} // namespace

TEST_CASE("confidence walk-up and minimum hit groups (khash and compact agree)") {
    const unsigned k = 31;
    std::mt19937_64 rng(23);
    std::vector<std::string> genomes;
    for(int g = 0; g < 4; ++g) genomes.push_back(random_dna(30000, rng));
    khash_t(p) *taxmap(kh_init(p));
    int khr;
    {
        khint_t ki = kh_put(p, taxmap, 1, &khr); kh_val(taxmap, ki) = 0;
        for(tax_t t = 2; t <= 5; ++t) {ki = kh_put(p, taxmap, t, &khr); kh_val(taxmap, ki) = 1;}
    }
    khash_t(c) *db(kh_init(c));
    {
        ClassifierGeneric<score::Lex> builder(db, spvec_t{}, k, k, 1, true, false, true, true);
        KmerKeys<score::Lex> enc(builder.enc_);
        for(int g = 0; g < 4; ++g)
            enc.for_each([&](u64 km) {
                const khint_t ki = kh_put(c, db, km, &khr);
                kh_val(db, ki) = (khr == 0 && kh_val(db, ki) != tax_t(2 + g)) ? 1: tax_t(2 + g);
            }, genomes[g].data(), genomes[g].size());
    }
    char tmpl[] = "/tmp/bns_conf_XXXXXX";
    const std::string fq = std::string(mkdtemp(tmpl)) + "/r.fq";
    {
        std::FILE *a = std::fopen(fq.data(), "w");
        const std::string chim = genomes[0].substr(5000, 100) + genomes[1].substr(9000, 50);
        std::fprintf(a, "@chimera\n%s\n+\n%s\n", chim.data(), std::string(150, 'I').data());
        std::fprintf(a, "@clean\n%s\n+\n%s\n", genomes[2].substr(1234, 150).data(), std::string(150, 'I').data());
        for(int r = 0; r < 200; ++r) {
            const std::string s = (r % 4 == 3) ? random_dna(150, rng): genomes[r % 4].substr(rng() % 29000, 150);
            std::fprintf(a, "@r%d\n%s\n+\n%s\n", r, s.data(), std::string(150, 'I').data());
        }
        std::fclose(a);
    }
    CompactKmerMap cm(db, CompactKmerMap::DEFAULT_MAX_LOAD, 2 * k);
    ClassifierGeneric<score::Lex> ck(db, spvec_t{}, k, k, 2, true, false, true, true);
    ClassifierGeneric<score::Lex, CompactKmerMap> cc(&cm, spvec_t{}, k, k, 2, true, false, true, true);
    auto both = [&](double conf, unsigned groups) {
        ck.set_confidence(conf); cc.set_confidence(conf);
        ck.set_min_hit_groups(groups); cc.set_min_hit_groups(groups);
        std::string a, b;
        auto ra = run_and_parse(ck, taxmap, fq, &a);
        auto rb = run_and_parse(cc, taxmap, fq, &b);
        REQUIRE(a == b);
        return ra;
    };
    const auto base = both(0., 0);
    REQUIRE(base.at("chimera") == 2);
    REQUIRE(base.at("clean") == 4);
    REQUIRE(both(0.5, 0).at("chimera") == 2);
    REQUIRE(both(0.7, 0).at("chimera") == 1);
    REQUIRE(both(0.8, 0).at("chimera") == 0);
    REQUIRE(both(0., 120).at("clean") == 4);
    REQUIRE(both(0., 121).at("clean") == 0);
    // confidence 1.0: every error-free read keeps its taxon, every random read is unclassified
    const auto full = both(1.0, 0);
    size_t kept = 0, rnd_called = 0;
    for(int r = 0; r < 200; ++r) {
        const std::string name = "r" + std::to_string(r);
        if(r % 4 == 3) rnd_called += full.at(name) != 0;
        else kept += full.at(name) == tax_t(2 + r % 4);
    }
    REQUIRE(kept == 150);
    REQUIRE(rnd_called == 0);
    // monotone: raising the threshold never classifies more reads
    size_t prev = ~size_t(0);
    for(const double conf: {0., 0.25, 0.5, 0.75, 0.9, 1.0}) {
        const auto r = both(conf, 0);
        size_t n = 0;
        for(const auto &kv: r) n += kv.second != 0;
        REQUIRE(n <= prev);
        prev = n;
    }
    REQUIRE_THROWS(ck.set_confidence(1.5));
    kh_destroy(c, db);
    kh_destroy(p, taxmap);
}
