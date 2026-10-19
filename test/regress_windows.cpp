// Windowed adjudication, lineage evidence and the second high-confidence call
// (classifier.h), on khash and compact databases, which must agree:
//  - taxonomy 1 -> families 10, 20; 10 -> species 2, 3; 20 -> species 4.
//    Species 4's genome carries a 1500 bp conserved island (rRNA/IS-like).
//  - a "novel relative" long read: 3500 bp of genome 2 at 10% substitutions
//    (few intact 31-mers, all hitting 2) + the island (exact hits to 4). The
//    pooled call is 4, a wrong-family call driven by the island alone; with
//    512-position windows the island wins only the windows it covers and the
//    read is called 2. E: shows the conflict against the pooled call;
//  - clean long reads keep their taxon with windows on; short reads (one
//    window) classify identically with and without windows;
//  - nearest-hit attribution (-N) also calls the novel read 2;
//  - H: (high-confidence call) is the call or an ancestor of it; with -X <= -x
//    it equals the call.
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
    u32 wcons = 0, wtot = 0;
    bool has_e = false, has_w = false, has_h = false;
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
            } else if(tok.rfind("W:", 0) == 0) {
                c.has_w = true;
                std::sscanf(tok.data() + 2, "%u/%u", &c.wcons, &c.wtot);
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

TEST_CASE("windowed adjudication rejects a conserved-island wrong-family call") {
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
    char tmpl[] = "/tmp/bns_win_XXXXXX";
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
    auto both = [&](const std::string &fq, unsigned window, double conf, double conf_hi, bool evidence) {
        for(auto *x: {&ck}) {x->set_window(window); x->set_confidence(conf); x->set_confidence_hi(conf_hi); x->set_emit_evidence(evidence);}
        cc.set_window(window); cc.set_confidence(conf); cc.set_confidence_hi(conf_hi); cc.set_emit_evidence(evidence);
        std::string a, b;
        auto ra = run_and_parse(ck, taxmap, fq, &a);
        auto rb = run_and_parse(cc, taxmap, fq, &b);
        REQUIRE(a == b);
        return ra;
    };

    // Pooled: the island drives a wrong-family call; E: shows the conflicting hits.
    const auto pooled = both(fq_long, 0, 0., 0., true);
    const Call &np = pooled.at("novel");
    std::fprintf(stderr, "[windows] novel read pooled: taxon %u E:%u,%u,%u,%u\n", np.taxon, np.e[0], np.e[1], np.e[2], np.e[3]);
    REQUIRE(np.taxon == 4);
    REQUIRE(np.has_e);
    REQUIRE(!np.has_w);
    REQUIRE(np.e[3] > 0);                // genome-2 hits conflict with the call
    REQUIRE(np.e[0] > np.e[3]);          // ... but are outnumbered k-mer for k-mer
    // Windowed: one vote per window -> called 2, and the W: field reports agreement.
    const auto windowed = both(fq_long, 512, 0., 0., true);
    const Call &nw = windowed.at("novel");
    std::fprintf(stderr, "[windows] novel read windowed: taxon %u W:%u/%u E:%u,%u,%u,%u\n",
                 nw.taxon, nw.wcons, nw.wtot, nw.e[0], nw.e[1], nw.e[2], nw.e[3]);
    REQUIRE(nw.taxon == 2);
    REQUIRE(nw.has_w);
    REQUIRE(nw.wtot == 9);               // 5000 bp -> 4970 positions -> 9 windows (tail merged)
    REQUIRE(2 * nw.wcons >= nw.wtot);
    REQUIRE(nw.wcons < nw.wtot);         // the island windows disagree
    // Nearest-hit attribution: the island owns only its own 1500 positions.
    ck.set_nearest(true); cc.set_nearest(true);
    {
        const auto r = both(fq_long, 0, 0., 0., false);
        REQUIRE(r.at("novel").taxon == 2);
        for(const tax_t t: {2, 3, 4})
            for(int i = 0; i < 5; ++i) REQUIRE(r.at("clean" + std::to_string(t) + "_" + std::to_string(i)).taxon == t);
        // error-free short reads from one genome: every hit agrees, so the same call as pooled
        const auto sh = both(fq_short, 0, 0., 0., false);
        ck.set_nearest(false); cc.set_nearest(false);
        const auto sp = both(fq_short, 0, 0., 0., false);
        for(int i = 0; i < 120; ++i)
            if(i % 4 != 1 && i % 10 != 7) REQUIRE(sh.at(short_names[i]).taxon == sp.at(short_names[i]).taxon);
    }
    ck.set_nearest(false); cc.set_nearest(false);
    // A strict agreement requirement (all windows) walks it up to the root.
    ck.set_window(512, 1.0); cc.set_window(512, 1.0);
    {
        std::string a, b;
        const auto ra = run_and_parse(ck, taxmap, fq_long, &a);
        const auto rb = run_and_parse(cc, taxmap, fq_long, &b);
        REQUIRE(a == b);
        REQUIRE(ra.at("novel").taxon == 1);
        for(const tax_t t: {2, 3, 4})
            for(int r = 0; r < 5; ++r) REQUIRE(ra.at("clean" + std::to_string(t) + "_" + std::to_string(r)).taxon == t);
    }
    // Clean long reads: unchanged by windows, and all evidence is exact.
    for(const tax_t t: {2, 3, 4})
        for(int r = 0; r < 5; ++r) {
            const std::string name = "clean" + std::to_string(t) + "_" + std::to_string(r);
            REQUIRE(pooled.at(name).taxon == t);
            REQUIRE(windowed.at(name).taxon == t);
            REQUIRE(windowed.at(name).wcons == windowed.at(name).wtot);
            REQUIRE(windowed.at(name).e[3] == 0);
        }
    // Short reads are one window: identical calls with and without windows.
    {
        const auto a = both(fq_short, 0, 0., 0., false), b = both(fq_short, 512, 0., 0., false);
        REQUIRE(a.size() == short_names.size());
        for(const auto &n: short_names) REQUIRE(a.at(n).taxon == b.at(n).taxon);
    }
    // Second call: ancestor-or-self of the call (or unclassified); equal to it when -X <= -x.
    {
        const auto r = both(fq_short, 0, 0.3, 0.9, true);
        size_t walked = 0;
        for(const auto &n: short_names) {
            const Call &c = r.at(n);
            REQUIRE(c.has_h);
            if(!c.taxon) {REQUIRE(c.hi == 0); continue;}
            REQUIRE((c.hi == 0 || ancestor_or_self(c.hi, c.taxon, taxmap)));
            walked += c.hi != c.taxon;
        }
        std::fprintf(stderr, "[windows] short reads: high-confidence call differs from the call on %zu/%zu\n", walked, short_names.size());
        REQUIRE(walked > 0);             // mutated reads fall short of 0.9
        const auto same = both(fq_short, 0, 0.5, 0.5, false);
        for(const auto &n: short_names) REQUIRE(same.at(n).hi == same.at(n).taxon);
    }
    REQUIRE_THROWS(ck.set_window(512, 1.5));
    REQUIRE_THROWS(ck.set_confidence_hi(-0.1));
    kh_destroy(c, db);
    kh_destroy(p, taxmap);
}
