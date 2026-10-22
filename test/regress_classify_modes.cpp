// Unit and end-to-end coverage of the classify modes in classifier.h / uniq.h:
//  - lineage_relation / lineage_evidence / confident_ancestor on a fixed tree,
//    with exact expectations (fractional totals included);
//  - identity-capped placement (-I) end to end: reads of a database strain at
//    1% / 8% / 16% substitutions stay at the strain / are capped at the family /
//    at the order; I: tracks 1 - divergence; custom thresholds; FASTA input
//    (no qualities: hits + misses) and paired-end reads give the same calls;
//    khash and compact databases agree;
//  - expected_intact_windows on lowercase and RNA bases;
//  - uniq report: pairs count once, unclassified reads land in row 0, a k = 40
//    exact 128-bit compact database and the khash database (61-bit hashed keys)
//    give identical exact counts, no coverage column without database counts;
//    UniqCounter merges across the exact/sketch boundary (and precisions).
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
std::string revcomp(const std::string &s) {
    std::string r(s.rbegin(), s.rend());
    for(auto &c: r) c = c == 'A' ? 'T': c == 'C' ? 'G': c == 'G' ? 'C': 'A';
    return r;
}
std::string slurp(std::FILE *f) {
    std::fflush(f); std::rewind(f);
    std::string out;
    char buf[1 << 14];
    size_t n;
    while((n = std::fread(buf, 1, sizeof(buf), f)) > 0) out.append(buf, n);
    return out;
}
struct Call {tax_t taxon = 0; double identity = -1.;};
std::map<std::string, Call> parse(const std::string &out) {
    std::map<std::string, Call> ret;
    std::istringstream is(out);
    std::string line;
    while(std::getline(is, line)) {
        std::istringstream ls(line);
        std::string cu, name, tok;
        Call c;
        ls >> cu >> name >> c.taxon;
        while(ls >> tok) if(tok.rfind("I:", 0) == 0) c.identity = std::stod(tok.substr(2));
        ret[name] = c;
    }
    return ret;
}
template<typename CL, typename Parents>
std::string run(const CL &cl, const Parents &taxmap, const std::string &fq1, const char *fq2 = nullptr) {
    std::FILE *o = std::tmpfile();
    process_dataset(cl, taxmap, fq1.data(), fq2, o, 1 << 16, 8);
    const std::string out = slurp(o);
    std::fclose(o);
    return out;
}
khash_t(p) *tree(const std::vector<std::pair<tax_t, tax_t>> &pc) {
    khash_t(p) *t(kh_init(p));
    int khr;
    for(const auto &e: pc) {const khint_t ki = kh_put(p, t, e.first, &khr); kh_val(t, ki) = e.second;}
    return t;
}
// 1 -> 2 (domain) -> 100 (order) -> 543 (family) -> 561 (genus) -> 562 (species) -> 9 (strain);
// 543 -> 620 -> 623; 100 -> 1903414 (family) -> 583 -> 584
const std::vector<std::pair<tax_t, tax_t>> TREE {
    {1, 0}, {2, 1}, {100, 2}, {543, 100}, {561, 543}, {562, 561}, {9, 562},
    {620, 543}, {623, 620}, {1903414, 100}, {583, 1903414}, {584, 583}};
rank_map_t tree_ranks() {
    return {{1, RANK_UNKNOWN}, {2, RANK_DOMAIN}, {100, RANK_ORDER}, {543, RANK_FAMILY}, {561, RANK_GENUS}, {562, RANK_SPECIES},
            {9, RANK_BELOW_SPECIES}, {620, RANK_GENUS}, {623, RANK_SPECIES}, {1903414, RANK_FAMILY}, {583, RANK_GENUS}, {584, RANK_SPECIES}};
}

} // namespace

TEST_CASE("lineage relation, evidence and confidence walk-up on a fixed tree") {
    khash_t(p) *kt = tree(TREE);
    const CompactTaxMap ct(kt);
    auto check = [](const auto &tm) {
    REQUIRE(lineage_relation(562, 562, tm) == 0);
    REQUIRE(lineage_relation(561, 562, tm) == 1);
    REQUIRE(lineage_relation(1, 562, tm) == 1);
    REQUIRE(lineage_relation(9, 562, tm) == 2);
    REQUIRE(lineage_relation(623, 562, tm) == 3);
    REQUIRE(lineage_relation(584, 9, tm) == 3);
    tax_counter hits;
    hits.add(562, 5); hits.add(561, 3); hits.add(9, 2); hits.add(623, 4);
    ReadEvidence ev;
    lineage_evidence(562, hits, tm, ev);
    REQUIRE((ev.exact == 5 && ev.ancestor == 3 && ev.descendant == 2 && ev.conflict == 4));
    lineage_evidence(0, hits, tm, ev);
    REQUIRE((ev.exact == 0 && ev.conflict == 14));
    // clades: 562 -> 7, 561 -> 10, 543 -> 14 (= everything)
    REQUIRE(confident_ancestor(9, hits, tm, 20., 0.1) == 9);      // 9's clade: 2 >= 2
    REQUIRE(confident_ancestor(9, hits, tm, 20., 0.35) == 562);   // 7 >= 7
    REQUIRE(confident_ancestor(9, hits, tm, 20., 0.5) == 561);    // 10 >= 10
    REQUIRE(confident_ancestor(9, hits, tm, 20., 0.7) == 543);    // 14 >= 14
    REQUIRE(confident_ancestor(9, hits, tm, 20., 0.75) == 0);     // even the root holds 14 < 15
    REQUIRE(confident_ancestor(9, hits, tm, 14., 1.0) == 543);    // hits-only denominator
    REQUIRE(confident_ancestor(9, hits, tm, 14.5, 0.96) == 543);  // fractional (quality-aware) total: 13.92 <= 14
    REQUIRE(confident_ancestor(9, hits, tm, 14.5, 0.99) == 0);    // 14.355 > 14
    };
    check(kt);   // khash_t(p)
    check(ct);   // CompactTaxMap: the same answers
    kh_destroy(p, kt);
}

TEST_CASE("expected_intact_windows: lowercase and RNA bases are valid") {
    const std::string upper(60, 'A'), lower(60, 'a'), rna(60, 'U'), q(60, char(33 + 30));
    const double e = expected_intact_windows(upper.data(), q.data(), 60, 31);
    REQUIRE(e == Approx(30 * std::pow(1. - 1e-3, 31)));
    REQUIRE(expected_intact_windows(lower.data(), q.data(), 60, 31) == Approx(e));
    REQUIRE(expected_intact_windows(rna.data(), q.data(), 60, 31) == Approx(e));
}

TEST_CASE("identity-capped placement end to end (khash and compact; FASTQ, FASTA, paired)") {
    const unsigned k = 31;
    std::mt19937_64 rng(77);
    khash_t(p) *tm = tree(TREE);
    const rank_map_t ranks = tree_ranks();
    const std::vector<std::pair<tax_t, std::string>> genomes {{9, random_dna(60000, rng)}, {623, random_dna(60000, rng)}, {584, random_dna(60000, rng)}};
    khash_t(c) *db(kh_init(c));
    {
        ClassifierGeneric<score::Lex> b(db, spvec_t{}, k, k, 1, true, false, true, true);
        KmerKeys<score::Lex> enc(b.enc_);
        int khr;
        for(const auto &g: genomes)
            enc.for_each([&](u64 km) {
                const khint_t ki = kh_put(c, db, km, &khr);
                kh_val(db, ki) = khr ? g.first: lca(tm, kh_val(db, ki), g.first);
            }, g.second.data(), g.second.size());
    }
    CompactKmerMap cm(db, CompactKmerMap::DEFAULT_MAX_LOAD, 2 * k);
    char tmpl[] = "/tmp/bns_modes_XXXXXX";
    const std::string dir(mkdtemp(tmpl));
    const std::string fq = dir + "/r.fq", fa = dir + "/r.fa", p1 = dir + "/r_1.fq", p2 = dir + "/r_2.fq";
    const double divs[] = {0.01, 0.08, 0.16};
    {
        std::FILE *a = std::fopen(fq.data(), "w"), *b = std::fopen(fa.data(), "w"), *m1 = std::fopen(p1.data(), "w"), *m2 = std::fopen(p2.data(), "w");
        for(int d = 0; d < 3; ++d)
            for(int r = 0; r < 4; ++r) {
                const std::string s = mutate(genomes[0].second.substr(rng() % 39000, 20000), divs[d], rng), q(s.size(), char(33 + 40));
                std::fprintf(a, "@d%d_%d\n%s\n+\n%s\n", d, r, s.data(), q.data());
                std::fprintf(b, ">d%d_%d\n%s\n", d, r, s.data());
                const std::string s1 = s.substr(0, 10000), s2 = revcomp(s.substr(10000));
                std::fprintf(m1, "@d%d_%d\n%s\n+\n%s\n", d, r, s1.data(), q.substr(0, 10000).data());
                std::fprintf(m2, "@d%d_%d\n%s\n+\n%s\n", d, r, s2.data(), q.substr(0, 10000).data());
            }
        std::fclose(a); std::fclose(b); std::fclose(m1); std::fclose(m2);
    }
    ClassifierGeneric<score::Lex> ck(db, spvec_t{}, k, k, 2, true, false, true, true);
    ClassifierGeneric<score::Lex, CompactKmerMap> cc(&cm, spvec_t{}, k, k, 2, true, false, true, true);
    // khash and compact databases, each with khash and CompactTaxMap parents: identical output
    const CompactTaxMap ctm(tm);
    auto both = [&](const std::string &in, const char *in2 = nullptr) {
        const std::string a = run(ck, tm, in, in2), b = run(cc, tm, in, in2);
        REQUIRE(a == b);
        REQUIRE(run(ck, ctm, in, in2) == a);
        REQUIRE(run(cc, ctm, in, in2) == a);
        return parse(a);
    };
    ck.set_emit_evidence(true); cc.set_emit_evidence(true);
    // without placement every read is called the strain (all hits are its k-mers)
    const auto plain = both(fq);
    // with placement: 1% stays at the strain, 8% is capped at the family, 16% at the order
    ck.set_identity_placement(&ranks); cc.set_identity_placement(&ranks);
    const auto placed = both(fq), fasta = both(fa), paired = both(p1, p2.data());
    const tax_t want[] = {9, 543, 100};
    for(int d = 0; d < 3; ++d)
        for(int r = 0; r < 4; ++r) {
            const std::string n = "d" + std::to_string(d) + "_" + std::to_string(r);
            std::fprintf(stderr, "[modes] %s: plain %u placed %u (I %.4f) fasta %u (I %.4f) paired %u (I %.4f)\n", n.data(),
                         plain.at(n).taxon, placed.at(n).taxon, placed.at(n).identity, fasta.at(n).taxon, fasta.at(n).identity,
                         paired.at(n).taxon, paired.at(n).identity);
            REQUIRE(plain.at(n).taxon == 9);
            REQUIRE(placed.at(n).taxon == want[d]);
            REQUIRE(std::abs(placed.at(n).identity - (1. - divs[d])) < 0.01);
            REQUIRE(fasta.at(n).taxon == want[d]);
            REQUIRE(std::abs(fasta.at(n).identity - placed.at(n).identity) < 0.002);
            REQUIRE(paired.at(n).taxon == want[d]);
            REQUIRE(std::abs(paired.at(n).identity - placed.at(n).identity) < 0.002);
        }
    // stricter thresholds: 1% divergence now falls below species, 8% below family
    ck.set_identity_placement(&ranks, 0.995, 0.95); cc.set_identity_placement(&ranks, 0.995, 0.95);
    const auto strict = both(fq);
    for(int r = 0; r < 4; ++r) {
        REQUIRE(strict.at("d0_" + std::to_string(r)).taxon == 543);
        REQUIRE(strict.at("d1_" + std::to_string(r)).taxon == 100);
    }
    // confidence walk-up (-x, -X) under both parent maps
    ck.set_identity_placement(nullptr); cc.set_identity_placement(nullptr);
    ck.set_confidence(0.6); cc.set_confidence(0.6); ck.set_confidence_hi(0.9); cc.set_confidence_hi(0.9);
    both(fq);
    REQUIRE_THROWS(ck.set_identity_placement(&ranks, 0.8, 0.9));
    REQUIRE_THROWS(ck.set_identity_placement(&ranks, 1.1, 0.9));
    kh_destroy(c, db);
    kh_destroy(p, tm);
}

TEST_CASE("uniq report: pairs, unclassified reads, 128-bit compact keys, no coverage") {
    std::mt19937_64 rng(5);
    khash_t(p) *tm = tree(TREE);
    const unsigned k = 40;   // khash: 61-bit hashed keys; compact 128: exact 80-bit keys
    const std::vector<std::pair<tax_t, std::string>> genomes {{562, random_dna(20000, rng)}, {584, random_dna(20000, rng)}};
    khash_t(c) *db(kh_init(c));
    ClassifierGeneric<score::Lex, CompactKmerMap128> b128(nullptr, spvec_t{}, k, k, 1, true, false, true, true);
    KmerKeys<score::Lex, u128> e128(b128.enc_);
    CompactLcaBuilder<u128> builder(80000, 16, KmerKeys<score::Lex, u128>::key_bits_for(k));
    {
        ClassifierGeneric<score::Lex> b(db, spvec_t{}, k, k, 1, true, false, true, true);
        KmerKeys<score::Lex> enc(b.enc_);
        int khr;
        for(const auto &g: genomes) {
            enc.for_each([&](u64 km) {const khint_t ki = kh_put(c, db, km, &khr); kh_val(db, ki) = g.first;}, g.second.data(), g.second.size());
            e128.for_each([&](u128 km) {builder.add(km, tm, g.first);}, g.second.data(), g.second.size());
        }
    }
    auto m128 = builder.finish();
    REQUIRE(m128->size() == kh_size(db));
    char tmpl[] = "/tmp/bns_uq_XXXXXX";
    const std::string dir(mkdtemp(tmpl)), p1 = dir + "/a_1.fq", p2 = dir + "/a_2.fq";
    {
        std::FILE *a = std::fopen(p1.data(), "w"), *c2 = std::fopen(p2.data(), "w");
        for(int r = 0; r < 2; ++r) {   // 2 pairs from 562 (mates 300 bp apart; 444 k-mers: exact counters), 3 random pairs
            const size_t off = rng() % 19000;
            const std::string m1 = genomes[0].second.substr(off, 150), m2 = revcomp(genomes[0].second.substr(off + 450, 150));
            std::fprintf(a, "@p%d\n%s\n+\n%s\n", r, m1.data(), std::string(150, 'I').data());
            std::fprintf(c2, "@p%d\n%s\n+\n%s\n", r, m2.data(), std::string(150, 'I').data());
        }
        for(int r = 0; r < 3; ++r) {
            const std::string m1 = random_dna(150, rng), m2 = random_dna(150, rng);
            std::fprintf(a, "@u%d\n%s\n+\n%s\n", r, m1.data(), std::string(150, 'I').data());
            std::fprintf(c2, "@u%d\n%s\n+\n%s\n", r, m2.data(), std::string(150, 'I').data());
        }
        std::fclose(a); std::fclose(c2);
    }
    auto report = [&](auto &cl) {
        cl.enable_uniq();
        run(cl, tm, p1, p2.data());
        std::FILE *u = std::tmpfile();
        write_uniq_report(u, *cl.uniq(), tm, nullptr, [](tax_t) {return std::string("-");});
        const std::string r = slurp(u);
        std::fclose(u);
        return r;
    };
    ClassifierGeneric<score::Lex> ck(db, spvec_t{}, k, k, 2, true, false, true, true);
    ClassifierGeneric<score::Lex, CompactKmerMap128> c8(m128.get(), spvec_t{}, k, k, 2, true, false, true, true);
    const std::string rk = report(ck), r8 = report(c8);
    std::fprintf(stderr, "[modes] paired uniq report:\n%s", rk.data());
    REQUIRE(rk == r8);   // exact counters (< 512 distinct): same counts from hashed and exact keys
    std::istringstream is(rk);
    std::string line;
    bool saw0 = false, saw562 = false;
    while(std::getline(is, line)) {
        if(line[0] == '#') continue;
        std::istringstream ls(line);
        tax_t t, par;
        std::string rank, cov;
        double pct, reads, treads, kmers, uniq, tuniq, dup;
        ls >> t >> par >> rank >> pct >> reads >> treads >> kmers >> uniq >> tuniq >> dup >> cov;
        if(t == 0) {saw0 = true; REQUIRE(reads == 3);}                     // 3 random pairs
        if(t == 562) {
            saw562 = true;
            REQUIRE(reads == 2);                                            // a pair counts once
            REQUIRE(kmers == 2 * 2 * (150 - k + 1));                       // both mates' k-mers
            REQUIRE(uniq <= kmers);
            REQUIRE(cov == "NA");                                           // no database counts given
        }
    }
    REQUIRE((saw0 && saw562));
    kh_destroy(c, db);
    kh_destroy(p, tm);
}

TEST_CASE("UniqCounter merges across the exact / sketch boundary") {
    UniqCounter a, b, c(10);
    for(uint64_t i = 0; i < 400; ++i) a.add(uniq_hash(i));
    for(uint64_t i = 200; i < 600; ++i) b.add(uniq_hash(i));
    REQUIRE((!a.sketched() && !b.sketched()));
    a.merge(b);                        // 600 distinct: past EXACT_MAX once compacted
    const double e = a.estimate();
    REQUIRE(std::abs(e - 600.) < 30.);
    REQUIRE(a.estimate() == e);        // estimate is stable
    for(uint64_t i = 0; i < 50000; ++i) c.add(uniq_hash(i + 7777777));
    REQUIRE(c.sketched());
    REQUIRE(std::abs(c.estimate() - 50000.) < 0.1 * 50000.);   // p = 10: ~3.2% standard error
    UniqCounter d;
    for(uint64_t i = 0; i < 10; ++i) d.add(uniq_hash(i + 7777777));
    d.merge(c);                        // exact (default p = 12) absorbs a p = 10 sketch: becomes a p = 10 sketch
    REQUIRE(d.sketched());
    REQUIRE(std::abs(d.estimate() - 50000.) < 0.1 * 50000.);
    UniqCounter e12;
    for(uint64_t i = 0; i < 5000; ++i) e12.add(uniq_hash(i));
    REQUIRE_THROWS(e12.merge(c));      // two sketches of different precision
}
