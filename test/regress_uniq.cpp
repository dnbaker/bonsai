// KrakenUniq-style per-taxon distinct k-mer report (classify -U; uniq.h):
//  - UniqCounter: exact below EXACT_MAX distinct hashes, HyperLogLog beyond
//    (estimate within 3% at 100k), merge = union;
//  - 60 distinct reads from genome 2 vs one read from genome 3 repeated 20
//    times: taxon 3 gets 20 reads but only 120 distinct k-mers (dup 20.0);
//    taxon 2's distinct count matches the exact union of its reads' k-mers
//    (within HyperLogLog error: it has > EXACT_MAX);
//  - coverage = distinct / database k-mers of the clade;
//  - 1 and 4 threads, khash and compact databases: identical reports;
//  - classification output is unchanged by -U.
#include "test/catch.hpp"

#include "classifier.h"
#include "kmermap.h"

#include <cstdio>
#include <map>
#include <random>
#include <set>
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
std::string slurp(std::FILE *f) {
    std::fflush(f); std::rewind(f);
    std::string out;
    char buf[1 << 14];
    size_t n;
    while((n = std::fread(buf, 1, sizeof(buf), f)) > 0) out.append(buf, n);
    return out;
}
struct Row {double reads, taxon_reads, kmers, uniq, taxon_uniq, dup; std::string cov;};
std::map<tax_t, Row> parse_report(const std::string &rep) {
    std::map<tax_t, Row> ret;
    std::istringstream is(rep);
    std::string line;
    while(std::getline(is, line)) {
        if(line.empty() || line[0] == '#') continue;
        std::istringstream ls(line);
        tax_t t, par;
        std::string rank;
        double pct;
        Row r;
        ls >> t >> par >> rank >> pct >> r.reads >> r.taxon_reads >> r.kmers >> r.uniq >> r.taxon_uniq >> r.dup >> r.cov;
        ret[t] = r;
    }
    return ret;
}

} // namespace

TEST_CASE("UniqCounter: exact, then HyperLogLog; merge is union") {
    UniqCounter a, b;
    for(uint64_t i = 0; i < 300; ++i) a.add(uniq_hash(i)), a.add(uniq_hash(i));
    REQUIRE(!a.sketched());
    REQUIRE(a.estimate() == 300.);
    for(uint64_t i = 0; i < 100000; ++i) b.add(uniq_hash(i + 1000000));
    REQUIRE(b.sketched());
    REQUIRE(std::abs(b.estimate() - 100000.) < 3000.);
    UniqCounter c;
    for(uint64_t i = 0; i < 200; ++i) c.add(uniq_hash(i + 150));   // 150 shared with a
    a.merge(c);
    REQUIRE(a.estimate() == 350.);
    a.merge(b);
    REQUIRE(a.sketched());
    REQUIRE(std::abs(a.estimate() - 100350.) < 3500.);
}

TEST_CASE("KrakenUniq-style report: distinct k-mers expose a repeated-read taxon") {
    const unsigned k = 31;
    std::mt19937_64 rng(41);
    std::vector<std::string> genomes {random_dna(30000, rng), random_dna(30000, rng)};
    khash_t(p) *taxmap(kh_init(p));
    int khr;
    for(const auto &pc: std::vector<std::pair<tax_t, tax_t>>{{1, 0}, {10, 1}, {2, 10}, {3, 10}}) {
        const khint_t ki = kh_put(p, taxmap, pc.first, &khr);
        kh_val(taxmap, ki) = pc.second;
    }
    khash_t(c) *db(kh_init(c));
    std::set<u64> keys2;   // exact distinct k-mers of genome-2 reads
    {
        ClassifierGeneric<score::Lex> builder(db, spvec_t{}, k, k, 1, true, false, true, true);
        KmerKeys<score::Lex> enc(builder.enc_);
        for(int g = 0; g < 2; ++g)
            enc.for_each([&](u64 km) {
                const khint_t ki = kh_put(c, db, km, &khr);
                kh_val(db, ki) = tax_t(2 + g);
            }, genomes[g].data(), genomes[g].size());
    }
    char tmpl[] = "/tmp/bns_uniq_XXXXXX";
    const std::string fq = std::string(mkdtemp(tmpl)) + "/r.fq";
    {
        ClassifierGeneric<score::Lex> builder(db, spvec_t{}, k, k, 1, true, false, true, true);
        KmerKeys<score::Lex> enc(builder.enc_);
        std::FILE *a = std::fopen(fq.data(), "w");
        for(int r = 0; r < 60; ++r) {
            const std::string s = genomes[0].substr(rng() % 29000, 150);
            enc.for_each([&](u64 km) {keys2.insert(km);}, s.data(), s.size());
            std::fprintf(a, "@a%d\n%s\n+\n%s\n", r, s.data(), std::string(150, 'I').data());
        }
        const std::string dup = genomes[1].substr(777, 150);
        for(int r = 0; r < 20; ++r) std::fprintf(a, "@b%d\n%s\n+\n%s\n", r, dup.data(), std::string(150, 'I').data());
        std::fclose(a);
    }
    std::unordered_map<tax_t, uint64_t> db_kmers;
    for(khiter_t ki = 0; ki != kh_end(db); ++ki) if(kh_exist(db, ki)) ++db_kmers[kh_val(db, ki)];
    auto rank_name = [](tax_t) {return std::string("-");};
    CompactKmerMap cm(db, CompactKmerMap::DEFAULT_MAX_LOAD, 2 * k);
    auto run = [&](auto &cl, bool uniq, std::string *report) {
        if(uniq) cl.enable_uniq();
        std::FILE *o = std::tmpfile();
        process_dataset(cl, taxmap, fq.data(), nullptr, o, 32, 8);
        const std::string out = slurp(o);
        std::fclose(o);
        if(uniq && report) {
            std::FILE *u = std::tmpfile();
            write_uniq_report(u, *cl.uniq(), taxmap, &db_kmers, rank_name);
            *report = slurp(u);
            std::fclose(u);
        }
        return out;
    };
    std::string r1, r4, rc, plain_out;
    {
        ClassifierGeneric<score::Lex> cl(db, spvec_t{}, k, k, 1, true, false, true, true);
        plain_out = run(cl, false, nullptr);
    }
    std::string out1, out4, outc;
    {ClassifierGeneric<score::Lex> cl(db, spvec_t{}, k, k, 1, true, false, true, true); out1 = run(cl, true, &r1);}
    {ClassifierGeneric<score::Lex> cl(db, spvec_t{}, k, k, 4, true, false, true, true); out4 = run(cl, true, &r4);}
    {ClassifierGeneric<score::Lex, CompactKmerMap> cl(&cm, spvec_t{}, k, k, 4, true, false, true, true); outc = run(cl, true, &rc);}
    REQUIRE(out1 == plain_out);          // -U leaves classification output unchanged
    REQUIRE(out4 == plain_out);
    REQUIRE(outc == plain_out);
    REQUIRE(r1 == r4);                   // thread merging
    REQUIRE(r1 == rc);                   // khash and compact
    std::fprintf(stderr, "[uniq] report:\n%s", r1.data());
    const auto rows = parse_report(r1);
    REQUIRE(rows.at(3).reads == 20);
    REQUIRE(rows.at(3).uniq == 120);
    REQUIRE(rows.at(3).kmers == 2400);
    REQUIRE(std::abs(rows.at(3).dup - 20.0) < 1e-9);
    REQUIRE(rows.at(2).reads == 60);
    // 6000+ distinct: past EXACT_MAX, so a p = 12 HyperLogLog estimate (~1.6% standard error)
    REQUIRE(std::abs(rows.at(2).uniq - double(keys2.size())) < 0.05 * double(keys2.size()));
    REQUIRE(rows.at(2).kmers == 60 * 120);
    REQUIRE(rows.at(10).reads == 80);
    REQUIRE(std::abs(rows.at(10).uniq - double(keys2.size() + 120)) < 0.05 * double(keys2.size()));
    REQUIRE(std::abs(std::stod(rows.at(3).cov) - 120. / double(db_kmers.at(3))) < 1e-6);
    kh_destroy(c, db);
    kh_destroy(p, taxmap);
}
