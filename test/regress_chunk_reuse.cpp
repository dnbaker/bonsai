// process_dataset reads records in chunks and reuses the record array across
// chunks (bseq_realloc_read). A chunk with more records than any earlier one
// used slots that bseq_read had allocated but never initialized (realloc
// garbage), so rekseq2bseq1 realloc'd a garbage name pointer: on real HiFi
// data (118- vs 120-read chunks) a few records per run came out with empty
// names or merged contents, and some runs aborted. It also never grew the
// array past bseq_read's capacity (4096 records).
// Here: 2 x 100 kb reads (chunk 1: 2 records), then 5000 x 40 bp (one chunk
// of 5000 records, past the 4096 capacity), then alternating 3 kb / 40 bp
// reads; every record must come out once, named, with its taxon, identically
// on every run and thread count.
#include "test/catch.hpp"

#include "classifier.h"

#include <cstdio>
#include <map>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using namespace bns;

TEST_CASE("chunked reading: later chunks with more records than earlier ones") {
    const unsigned k = 31;
    std::mt19937_64 rng(5);
    std::string g(300000, 'A');
    for(auto &c: g) c = "ACGT"[rng() & 3];
    khash_t(p) *taxmap(kh_init(p));
    int khr;
    {khint_t ki = kh_put(p, taxmap, 1, &khr); kh_val(taxmap, ki) = 0; ki = kh_put(p, taxmap, 2, &khr); kh_val(taxmap, ki) = 1;}
    khash_t(c) *db(kh_init(c));
    {
        ClassifierGeneric<score::Lex> builder(db, spvec_t{}, k, k, 1, true, false, true, true);
        KmerKeys<score::Lex> enc(builder.enc_);
        enc.for_each([&](u64 km) {const khint_t ki = kh_put(c, db, km, &khr); kh_val(db, ki) = 2;}, g.data(), g.size());
    }
    char tmpl[] = "/tmp/bns_chunk_XXXXXX";
    const std::string dir(mkdtemp(tmpl)), fq = dir + "/r.fq";
    std::vector<std::string> names;
    {
        std::FILE *f = std::fopen(fq.data(), "w");
        auto put = [&](size_t len) {
            const size_t o = rng() % (g.size() - len);
            names.push_back("read" + std::to_string(names.size()));
            std::fprintf(f, "@%s extra/comment\n%s\n+\n%s\n", names.back().data(), g.substr(o, len).data(), std::string(len, 'I').data());
        };
        put(100000); put(100000);
        for(int i = 0; i < 5000; ++i) put(40);
        for(int i = 0; i < 400; ++i) put(i & 1 ? 40: 3000);
        std::fclose(f);
    }
    std::string first;
    for(const unsigned nt: {1u, 4u, 1u, 4u}) {
        ClassifierGeneric<score::Lex> cl(db, spvec_t{}, k, k, nt, true, false, true, true);
        std::FILE *o = std::tmpfile();
        process_dataset(cl, taxmap, fq.data(), nullptr, o, 200000, 32);
        std::fflush(o); std::rewind(o);
        std::string out;
        char buf[1 << 14];
        for(size_t n; (n = std::fread(buf, 1, sizeof(buf), o)) > 0;) out.append(buf, n);
        std::fclose(o);
        std::map<std::string, int> seen;
        std::istringstream is(out);
        std::string line;
        size_t lines = 0, bad = 0;
        while(std::getline(is, line)) {
            ++lines;
            std::istringstream ls(line);
            std::string cu, name; tax_t t = 0;
            ls >> cu >> name >> t;
            ++seen[name];
            bad += t != 2;
        }
        std::fprintf(stderr, "[chunk] threads %u: %zu lines, %zu distinct names, %zu not taxon 2\n", nt, lines, seen.size(), bad);
        REQUIRE(lines == names.size());
        REQUIRE(seen.size() == names.size());
        for(const auto &n: names) REQUIRE(seen[n] == 1);
        REQUIRE(bad == 0);
        if(first.empty()) first = out;
        else REQUIRE(out == first);
    }
    std::remove(fq.data()); rmdir(dir.data());
    kh_destroy(c, db);
    kh_destroy(p, taxmap);
}
