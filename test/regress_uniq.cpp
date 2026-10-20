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
#include <thread>
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

TEST_CASE("b-bit MinHash counters (-B): estimates, union, written sketches, containment") {
    UniqCounter a(10, UniqSketch::BBMH), b(10, UniqSketch::BBMH);
    for(uint64_t i = 0; i < 300; ++i) a.add(uniq_hash(i));
    REQUIRE(!a.sketched());
    REQUIRE(a.estimate() == 300.);
    for(uint64_t i = 0; i < 60000; ++i) b.add(uniq_hash(i));
    REQUIRE(b.sketched());
    REQUIRE(b.kind() == UniqSketch::BBMH);
    REQUIRE(std::abs(b.estimate() - 60000.) < 0.08 * 60000.);
    a.merge(b);                                    // subset absorbed: still ~60000
    REQUIRE(a.sketched());
    REQUIRE(std::abs(a.estimate() - 60000.) < 0.08 * 60000.);
    UniqCounter h(10);                             // a HyperLogLog sketch cannot merge into a MinHash one
    for(uint64_t i = 0; i < 5000; ++i) h.add(uniq_hash(i));
    REQUIRE_THROWS(a.merge(h));
    // written sketches: a taxon's read k-mers are contained in its genome's k-mers
    char tmpl[] = "/tmp/bns_bbmh_XXXXXX";
    const std::string dir(mkdtemp(tmpl));
    UniqCounter reads(10, UniqSketch::BBMH), genome(10, UniqSketch::BBMH), other(10, UniqSketch::BBMH);
    for(uint64_t i = 0; i < 200000; ++i) genome.add(uniq_hash(i));
    for(uint64_t i = 0; i < 50000; ++i) reads.add(uniq_hash(i * 4));        // a quarter of the genome
    for(uint64_t i = 0; i < 50000; ++i) other.add(uniq_hash(i + 10000000)); // disjoint
    reads.bbmh(16).write((dir + "/reads.bbmh").data());
    genome.bbmh(16).write((dir + "/genome.bbmh").data());
    const sketch::FinalBBitMinHash fr((dir + "/reads.bbmh").data()), fg((dir + "/genome.bbmh").data()), fo(other.bbmh(16));
    REQUIRE(std::abs(fr.cardinality_estimate() - 50000.) < 0.08 * 50000.);
    const double cont = fr.containment_index(fg), cont_o = fo.containment_index(fg);
    std::fprintf(stderr, "[uniq] bbmh containment: reads in genome %.3f, disjoint set in genome %.3f, Jaccard %.3f\n",
                 cont, cont_o, fr.jaccard_index(fg));
    REQUIRE(cont > 0.8);                           // p = 10: ~1024 registers, a few % noise
    REQUIRE(cont_o < 0.05);
    REQUIRE(std::abs(fr.jaccard_index(fg) - 0.25) < 0.05);
    // a sparse sketch (empty registers) keeps an accurate cardinality when finalized and written
    UniqCounter sparse(10, UniqSketch::BBMH);
    for(uint64_t i = 0; i < 3000; ++i) sparse.add(uniq_hash(i + 555));
    REQUIRE(sparse.sketched());                     // sparse: ~5% of the 1024 registers still empty
    REQUIRE(std::abs(sparse.estimate() - 3000.) < 0.1 * 3000.);
    sparse.bbmh(16).write((dir + "/sparse.bbmh").data());
    const sketch::FinalBBitMinHash fs((dir + "/sparse.bbmh").data());
    REQUIRE(std::abs(fs.cardinality_estimate() - 3000.) < 0.1 * 3000.);
    // and in the small range, where linear counting applies
    UniqCounter small(10, UniqSketch::BBMH);
    for(uint64_t i = 0; i < 1500; ++i) small.add(uniq_hash(i + 999));   // sketched at the 1024-hash compaction
    REQUIRE(small.sketched());
    REQUIRE(std::abs(small.estimate() - 1500.) < 0.1 * 1500.);
    REQUIRE(std::abs(small.bbmh(16).cardinality_estimate() - 1500.) < 0.1 * 1500.);
}

TEST_CASE("-B through the classifier: per-clade sketches written next to the report") {
    const unsigned k = 31;
    std::mt19937_64 rng(43);
    std::vector<std::string> genomes {random_dna(30000, rng), random_dna(30000, rng)};
    khash_t(p) *taxmap(kh_init(p));
    int khr;
    for(const auto &pc: std::vector<std::pair<tax_t, tax_t>>{{1, 0}, {10, 1}, {2, 10}, {3, 10}}) {
        const khint_t ki = kh_put(p, taxmap, pc.first, &khr);
        kh_val(taxmap, ki) = pc.second;
    }
    khash_t(c) *db(kh_init(c));
    {
        ClassifierGeneric<score::Lex> builder(db, spvec_t{}, k, k, 1, true, false, true, true);
        KmerKeys<score::Lex> enc(builder.enc_);
        for(int g = 0; g < 2; ++g)
            enc.for_each([&](u64 km) {const khint_t ki = kh_put(c, db, km, &khr); kh_val(db, ki) = tax_t(2 + g);}, genomes[g].data(), genomes[g].size());
    }
    char tmpl[] = "/tmp/bns_bbmhc_XXXXXX";
    const std::string dir(mkdtemp(tmpl)), fq = dir + "/r.fq";
    {
        std::FILE *a = std::fopen(fq.data(), "w");
        for(int r = 0; r < 100; ++r) {
            const std::string s = genomes[r % 2].substr(rng() % 29000, 150);
            std::fprintf(a, "@r%d\n%s\n+\n%s\n", r, s.data(), std::string(150, 'I').data());
        }
        std::fclose(a);
    }
    ClassifierGeneric<score::Lex> hl(db, spvec_t{}, k, k, 2, true, false, true, true), bm(db, spvec_t{}, k, k, 2, true, false, true, true);
    hl.enable_uniq(12, UniqSketch::HLL);
    bm.enable_uniq(10, UniqSketch::BBMH);
    for(auto *cl: {&hl, &bm}) {
        std::FILE *o = std::tmpfile();
        process_dataset(*cl, taxmap, fq.data(), nullptr, o, 32, 8);
        std::fclose(o);
    }
    auto rank_name = [](tax_t) {return std::string("-");};
    std::FILE *u1 = std::tmpfile(), *u2 = std::tmpfile();
    write_uniq_report(u1, *hl.uniq(), taxmap, nullptr, rank_name);
    write_uniq_report(u2, *bm.uniq(), taxmap, nullptr, rank_name, &dir, 16);
    const auto r1 = parse_report(slurp(u1)), r2 = parse_report(slurp(u2));
    std::fclose(u1); std::fclose(u2);
    for(const tax_t t: {1u, 10u, 2u, 3u}) {
        REQUIRE(r1.at(t).reads == r2.at(t).reads);
        REQUIRE(r1.at(t).kmers == r2.at(t).kmers);
        REQUIRE(std::abs(r1.at(t).uniq - r2.at(t).uniq) < 0.1 * r1.at(t).uniq);   // HLL vs MinHash estimates
        const sketch::FinalBBitMinHash f((dir + "/" + std::to_string(t) + ".bbmh").data());
        REQUIRE(std::abs(f.cardinality_estimate() - r2.at(t).uniq) < 0.1 * r2.at(t).uniq + 1);
    }
    kh_destroy(c, db);
    kh_destroy(p, taxmap);
}

TEST_CASE("UniqTable: concurrent batches give the serial counts (exact tallies locked, sketches lock-free)") {
    // Taxon t gets hashes [0, n_t) from every thread, in a different order per
    // thread; taxa 1-3 stay exact (<= 512 distinct), 4-6 become sketches part
    // way through, while other threads may still be in their exact phase.
    const std::vector<uint64_t> n{0, 50, 300, 512, 900, 5000, 60000};
    for(const UniqSketch kind: {UniqSketch::HLL, UniqSketch::BBMH}) {
        const unsigned p = kind == UniqSketch::HLL ? 12: 10, nt = 8;
        UniqTable shared(nt, p, kind), serial(1, p, kind);
        std::vector<std::thread> threads;
        for(unsigned tid = 0; tid < nt; ++tid)
            threads.emplace_back([&, tid] {
                std::mt19937_64 rng(tid);
                std::vector<std::pair<tax_t, uint64_t>> adds;
                for(tax_t t = 1; t < n.size(); ++t) for(uint64_t i = 0; i < n[t]; ++i) adds.emplace_back(t, uniq_hash(i * 7 + t));
                std::shuffle(adds.begin(), adds.end(), rng);
                UniqBatch &b = shared.batch(tid);
                for(const auto &a: adds) b.add_kmer(a.first, a.second);
                for(tax_t t = 1; t < n.size(); ++t) b.add_read(t);
            });
        for(auto &th: threads) th.join();
        for(unsigned tid = 0; tid < nt; ++tid)
            for(tax_t t = 1; t < n.size(); ++t) {
                for(uint64_t i = 0; i < n[t]; ++i) serial.batch(0).add_kmer(t, uniq_hash(i * 7 + t));
                serial.batch(0).add_read(t);
            }
        shared.flush(); serial.flush();
        for(tax_t t = 1; t < n.size(); ++t) {
            SharedTally &a = shared.at(t), &b = serial.at(t);
            REQUIRE(a.reads.load() == nt);
            REQUIRE(a.kmers.load() == nt * n[t]);
            REQUIRE(b.kmers.load() == a.kmers.load());
            REQUIRE(a.uniq.sketched() == (n[t] > UniqCounter::EXACT_MAX));
            const double ea = a.uniq.estimate(), eb = b.uniq.estimate();
            std::fprintf(stderr, "[uniq] table %s taxon %u: n=%llu shared %.1f serial %.1f\n", kind == UniqSketch::HLL ? "hll": "bbmh",
                         t, (unsigned long long)n[t], ea, eb);
            REQUIRE(ea == eb);   // register max/min are order-independent: identical sketches
            if(!a.uniq.sketched()) REQUIRE(ea == double(n[t]));
            else REQUIRE(std::abs(ea - double(n[t])) < 0.1 * double(n[t]));
        }
    }
}
