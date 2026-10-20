// Lineage evidence (-E), the second high-confidence call (-X) and the
// confidence denominators (default hits + misses, -R hits only, -Q expected
// error-free k-mers from base qualities), on khash and compact databases,
// which must agree:
//  - taxonomy 1 -> families 10, 20; 10 -> species 2, 3; 20 -> species 4.
//    Species 4's genome carries a 1500 bp conserved island (rRNA/IS-like).
//  - a "novel relative" long read: 3500 bp of genome 2 at 10% substitutions
//    (few intact 31-mers, all hitting 2) + the island (exact hits to 4). The
//    pooled call is 4, a wrong-family call; E: shows the conflicting hits;
//  - noisy reads from a database genome (5% substitutions, Q13): Kraken2's
//    confidence (hits + misses) leaves them unclassified, -Q keeps them, and
//    -Q still walks the novel read out of the wrong family (its hits fall far
//    short of the k-mers its qualities predict), while -R keeps the wrong call;
//  - H: is the call or an ancestor of it (or unclassified); with -X <= -x it
//    equals the call.
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
std::string mutate(std::string s, double rate, std::mt19937_64 &rng) {
    static const char b[] = "ACGT";
    std::uniform_real_distribution<double> u;
    for(auto &c: s) if(u(rng) < rate) {char n; do n = b[rng() & 3]; while(n == c); c = n;}
    return s;
}
struct Call {
    tax_t taxon = 0, hi = 0;
    u32 e[4] = {0, 0, 0, 0};
    bool has_e = false, has_h = false;
};
template<typename CL>
std::map<std::string, Call> run_and_parse(const CL &cl, const khash_t(p) *taxmap, const std::string &fq, std::string *raw = nullptr) {
    std::FILE *o = std::tmpfile();
    process_dataset(cl, taxmap, fq.data(), nullptr, o, 256, 32);
    std::fflush(o); std::rewind(o);
    std::string out;
    char buf[1 << 14];
    size_t n;
    while((n = std::fread(buf, 1, sizeof(buf), o)) > 0) out.append(buf, n);
    std::fclose(o);
    if(raw) *raw = out;
    std::map<std::string, Call> ret;
    std::istringstream is(out);
    std::string line;
    while(std::getline(is, line)) {
        std::istringstream ls(line);
        std::string cu, name, tok;
        Call c;
        ls >> cu >> name >> c.taxon;
        while(ls >> tok) {
            if(tok.rfind("E:", 0) == 0) {
                c.has_e = true;
                std::sscanf(tok.data() + 2, "%u,%u,%u,%u", &c.e[0], &c.e[1], &c.e[2], &c.e[3]);
            } else if(tok.rfind("H:", 0) == 0) {
                c.has_h = true;
                c.hi = tax_t(std::strtoul(tok.data() + 2, nullptr, 10));
            }
        }
        ret[name] = c;
    }
    return ret;
}
bool ancestor_or_self(tax_t a, tax_t t, const khash_t(p) *taxmap) {
    for(unsigned guard = 0; t && guard < 64; ++guard) {
        if(t == a) return true;
        t = get_parent(taxmap, t);
    }
    return false;
}

} // namespace

TEST_CASE("lineage evidence, second call and confidence denominators") {
    const unsigned k = 31;
    std::mt19937_64 rng(31);
    const std::string island = random_dna(1500, rng);
    std::map<tax_t, std::string> genomes;
    genomes[2] = random_dna(30000, rng);
    genomes[3] = random_dna(30000, rng);
    genomes[4] = random_dna(15000, rng) + island + random_dna(15000, rng);
    khash_t(p) *taxmap(kh_init(p));
    int khr;
    for(const auto &pc: std::vector<std::pair<tax_t, tax_t>>{{1, 0}, {10, 1}, {20, 1}, {2, 10}, {3, 10}, {4, 20}}) {
        const khint_t ki = kh_put(p, taxmap, pc.first, &khr);
        kh_val(taxmap, ki) = pc.second;
    }
    khash_t(c) *db(kh_init(c));
    {
        ClassifierGeneric<score::Lex> builder(db, spvec_t{}, k, k, 1, true, false, true, true);
        KmerKeys<score::Lex> enc(builder.enc_);
        for(const auto &g: genomes)
            enc.for_each([&](u64 km) {
                const khint_t ki = kh_put(c, db, km, &khr);
                kh_val(db, ki) = (khr == 0 && kh_val(db, ki) != g.first) ? 1: g.first;
            }, g.second.data(), g.second.size());
    }
    char tmpl[] = "/tmp/bns_evid_XXXXXX";
    const std::string dir(mkdtemp(tmpl));
    const std::string fq_long = dir + "/long.fq", fq_short = dir + "/short.fq";
    std::vector<std::string> short_names;
    {
        std::FILE *a = std::fopen(fq_long.data(), "w");
        const std::string novel = mutate(genomes[2].substr(2000, 3500), 0.10, rng) + island;
        std::fprintf(a, "@novel\n%s\n+\n%s\n", novel.data(), std::string(novel.size(), 'I').data());
        for(const tax_t t: {2, 3, 4})
            for(int r = 0; r < 5; ++r) {
                const std::string s = genomes[t].substr(rng() % 24000, 5000);
                std::fprintf(a, "@clean%u_%d\n%s\n+\n%s\n", t, r, s.data(), std::string(s.size(), 'I').data());
            }
        std::fclose(a);
        a = std::fopen(fq_short.data(), "w");
        for(int r = 0; r < 120; ++r) {
            const tax_t t = 2 + r % 3;
            std::string s = genomes[t].substr(rng() % 29000, 150);
            if(r % 4 == 1) s = mutate(s, 0.03, rng);
            if(r % 10 == 7) s = random_dna(150, rng);
            short_names.push_back("s" + std::to_string(r));
            std::fprintf(a, "@s%d\n%s\n+\n%s\n", r, s.data(), std::string(150, 'I').data());
        }
        std::fclose(a);
    }
    CompactKmerMap cm(db, CompactKmerMap::DEFAULT_MAX_LOAD, 2 * k);
    ClassifierGeneric<score::Lex> ck(db, spvec_t{}, k, k, 2, true, false, true, true);
    ClassifierGeneric<score::Lex, CompactKmerMap> cc(&cm, spvec_t{}, k, k, 2, true, false, true, true);
    auto both = [&](const std::string &fq, double conf, double conf_hi, bool evidence, bool hits_only = false, bool quality = false) {
        ck.set_confidence(conf); ck.set_confidence_hi(conf_hi); ck.set_emit_evidence(evidence);
        ck.set_confidence_hits_only(hits_only); ck.set_confidence_quality(quality);
        cc.set_confidence(conf); cc.set_confidence_hi(conf_hi); cc.set_emit_evidence(evidence);
        cc.set_confidence_hits_only(hits_only); cc.set_confidence_quality(quality);
        std::string a, b;
        auto ra = run_and_parse(ck, taxmap, fq, &a);
        auto rb = run_and_parse(cc, taxmap, fq, &b);
        REQUIRE(a == b);
        return ra;
    };
    // Pooled: the island drives a wrong-family call; E: shows the conflicting hits.
    const auto pooled = both(fq_long, 0., 0., true);
    const Call &np = pooled.at("novel");
    std::fprintf(stderr, "[evidence] novel read: taxon %u E:%u,%u,%u,%u\n", np.taxon, np.e[0], np.e[1], np.e[2], np.e[3]);
    REQUIRE(np.taxon == 4);
    REQUIRE(np.has_e);
    REQUIRE(np.e[3] > 0);                // genome-2 hits conflict with the call
    REQUIRE(np.e[0] > np.e[3]);          // ... but are outnumbered k-mer for k-mer
    // Clean long reads: all evidence is exact.
    for(const tax_t t: {2, 3, 4})
        for(int r = 0; r < 5; ++r) {
            const Call &c = pooled.at("clean" + std::to_string(t) + "_" + std::to_string(r));
            REQUIRE(c.taxon == t);
            REQUIRE(c.e[3] == 0);
            REQUIRE(c.e[0] + c.e[1] + c.e[2] > 0);
        }
    // Confidence denominators on noisy reads (5% substitutions, honest Q13).
    {
        const std::string fq_noisy = dir + "/noisy.fq";
        std::FILE *a = std::fopen(fq_noisy.data(), "w");
        const std::string novel = mutate(genomes[2].substr(2000, 3500), 0.10, rng) + island;
        std::fprintf(a, "@novel\n%s\n+\n%s\n", novel.data(), std::string(novel.size(), char(33 + 40)).data());
        for(const tax_t t: {2, 3, 4})
            for(int r = 0; r < 5; ++r) {
                const std::string s = mutate(genomes[t].substr(rng() % 24000, 5000), 0.05, rng);
                std::fprintf(a, "@noisy%u_%d\n%s\n+\n%s\n", t, r, s.data(), std::string(s.size(), char(33 + 13)).data());
            }
        std::fclose(a);
        const auto kraken = both(fq_noisy, 0.3, 0., false), hits = both(fq_noisy, 0.3, 0., false, true),
                   qual = both(fq_noisy, 0.3, 0., false, false, true);
        size_t k_kept = 0, h_kept = 0, q_kept = 0;
        for(const tax_t t: {2, 3, 4})
            for(int r = 0; r < 5; ++r) {
                const std::string n = "noisy" + std::to_string(t) + "_" + std::to_string(r);
                k_kept += kraken.at(n).taxon == t; h_kept += hits.at(n).taxon == t; q_kept += qual.at(n).taxon == t;
            }
        std::fprintf(stderr, "[evidence] noisy reads kept at 0.3: kraken %zu/15, hits-only %zu/15, quality %zu/15; novel read: %u / %u / %u\n",
                     k_kept, h_kept, q_kept, kraken.at("novel").taxon, hits.at("novel").taxon, qual.at("novel").taxon);
        REQUIRE(k_kept == 0);                // ~0.95^31 = 20% of k-mers survive: below 0.3 of hits + misses
        REQUIRE(h_kept == 15);
        REQUIRE(q_kept == 15);
        REQUIRE(hits.at("novel").taxon == 4);   // hits-only keeps the island call (90% of hits)
        // quality-aware: 1636 hits of ~4970 expected at Q40; 4's clade (1470) falls short of 0.3,
        // the root (0.33) does not: walked up out of the wrong family, as Kraken2 does on this read
        REQUIRE(qual.at("novel").taxon == 1);
        REQUIRE(kraken.at("novel").taxon == 1);
    }
    // Second call: ancestor-or-self of the call (or unclassified); equal to it when -X <= -x.
    {
        const auto r = both(fq_short, 0.3, 0.9, true);
        size_t walked = 0;
        for(const auto &n: short_names) {
            const Call &c = r.at(n);
            REQUIRE(c.has_h);
            if(!c.taxon) {REQUIRE(c.hi == 0); continue;}
            REQUIRE((c.hi == 0 || ancestor_or_self(c.hi, c.taxon, taxmap)));
            walked += c.hi != c.taxon;
        }
        std::fprintf(stderr, "[evidence] short reads: high-confidence call differs from the call on %zu/%zu\n", walked, short_names.size());
        REQUIRE(walked > 0);             // mutated reads fall short of 0.9
        const auto same = both(fq_short, 0.5, 0.5, false);
        for(const auto &n: short_names) REQUIRE(same.at(n).hi == same.at(n).taxon);
    }
    REQUIRE_THROWS(ck.set_confidence_hi(-0.1));
    kh_destroy(c, db);
    kh_destroy(p, taxmap);
}

TEST_CASE("expected error-free windows from base qualities") {
    const std::string seq(100, 'A');
    const unsigned span = 31;
    // constant Q20 (p = 0.01): 70 windows, each intact with probability 0.99^31
    const std::string q20(100, char(33 + 20));
    REQUIRE(std::abs(expected_intact_windows(seq.data(), q20.data(), 100, span) - 70 * std::pow(0.99, 31)) < 1e-9);
    // an N at position 50 removes the 31 windows covering it
    std::string withn = seq;
    withn[50] = 'N';
    REQUIRE(std::abs(expected_intact_windows(withn.data(), q20.data(), 100, span) - 39 * std::pow(0.99, 31)) < 1e-9);
    // quality 0 counts as certainly wrong; Q93 is ~error-free
    std::string mixed(100, char(33 + 93));
    mixed[0] = char(33);
    REQUIRE(std::abs(expected_intact_windows(seq.data(), mixed.data(), 100, span) - 69.) < 1e-5);
    REQUIRE(expected_intact_windows(seq.data(), q20.data(), 30, span) == 0.);
}
