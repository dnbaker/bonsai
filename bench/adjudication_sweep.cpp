// Read-level adjudication modes on the conserved-island failure (Perseus,
// Nguyen & Schatz 2026: rRNA/IS regions shared across families drive
// wrong-family calls), on synthetic data:
//   taxonomy 1 -> families 10, 20; 10 -> species 2, 3; 20 -> species 4.
//   Species 4's genome carries an island of I bp. A read from an unsampled
//   relative of 2 is (L - I) bp of genome 2 at divergence d, with the island
//   spliced in at a random point. Truth: 2 (in-lineage: 10).
// Modes: pooled (default), Kraken2 confidence 0.3, hits-only confidence (-R)
// 0.5/0.7/0.9, quality-aware confidence (-Q) 0.3/0.5/0.7 (reads at Q40).
// (Windowed adjudication and nearest-hit attribution were tried here and on
// the real-genome benchmarks and dropped: no gain on real data.) For each
// (I, d, L) cell: reads // called 2, in lineage (2 or 10), wrong family (4 or 20), root/unclassified.
// Then single-thread throughput of each mode on short and long error-free reads.
// Usage: adjudication_sweep [reads_per_cell=20] [seed=7] > sweep.csv
#include "classifier.h"
#include "kmermap.h"

#include <chrono>
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
struct Mode {const char *name; double conf; bool hits_only, quality;};
const Mode MODES[] = {
    {"pooled", 0., false, false},
    {"conf0.3", 0.3, false, false},
    {"hits0.5", 0.5, true, false},
    {"hits0.7", 0.7, true, false},
    {"hits0.9", 0.9, true, false},
    {"qual0.3", 0.3, false, true},
    {"qual0.5", 0.5, false, true},
    {"qual0.7", 0.7, false, true},
};
template<typename CL>
void set_mode(CL &c, const Mode &m) {
    c.set_confidence(m.conf); c.set_confidence_hits_only(m.hits_only); c.set_confidence_quality(m.quality);
}
template<typename CL>
std::vector<tax_t> classify(const CL &cl, const khash_t(p) *taxmap, const std::string &fq) {
    std::FILE *o = std::tmpfile();
    process_dataset(cl, taxmap, fq.data(), nullptr, o, 1 << 14, 32);
    std::fflush(o); std::rewind(o);
    std::string out;
    char buf[1 << 16];
    size_t n;
    while((n = std::fread(buf, 1, sizeof(buf), o)) > 0) out.append(buf, n);
    std::fclose(o);
    std::vector<tax_t> ret;
    std::istringstream is(out);
    std::string line, cu, name;
    while(std::getline(is, line)) {
        std::istringstream ls(line);
        tax_t t = 0;
        ls >> cu >> name >> t;
        ret.push_back(t);
    }
    return ret;
}
void write_fq(const std::string &path, const std::vector<std::string> &reads) {
    std::FILE *a = std::fopen(path.data(), "w");
    for(size_t i = 0; i < reads.size(); ++i)
        std::fprintf(a, "@r%zu\n%s\n+\n%s\n", i, reads[i].data(), std::string(reads[i].size(), 'I').data());
    std::fclose(a);
}

} // namespace

int main(int argc, char **argv) {
    const unsigned per_cell = argc > 1 ? std::atoi(argv[1]): 20;
    const uint64_t seed = argc > 2 ? std::strtoull(argv[2], nullptr, 10): 7;
    const unsigned k = 31;
    std::fprintf(stderr, "[sweep] seed %llu, %u reads per cell, k = %u\n", (unsigned long long)seed, per_cell, k);
    std::mt19937_64 rng(seed);
    khash_t(p) *taxmap(kh_init(p));
    int khr;
    for(const auto &pc: std::vector<std::pair<tax_t, tax_t>>{{1, 0}, {10, 1}, {20, 1}, {2, 10}, {3, 10}, {4, 20}}) {
        const khint_t ki = kh_put(p, taxmap, pc.first, &khr);
        kh_val(taxmap, ki) = pc.second;
    }
    const std::string g2 = random_dna(40000, rng), g3 = random_dna(40000, rng);
    const std::string g4a = random_dna(20000, rng), g4b = random_dna(20000, rng);
    const std::string island_full = random_dna(3000, rng);
    char tmpl[] = "/tmp/bns_sweep_XXXXXX";
    const std::string dir(mkdtemp(tmpl)), fq = dir + "/r.fq";
    std::printf("island_bp,divergence,read_bp,mode,reads,called_true,in_lineage,wrong_family,root_or_unclassified\n");
    for(const unsigned isz: {0u, 300u, 1000u, 3000u}) {
        const std::string island = island_full.substr(0, isz);
        khash_t(c) *db(kh_init(c));
        {
            ClassifierGeneric<score::Lex> builder(db, spvec_t{}, k, k, 1, true, false, true, true);
            KmerKeys<score::Lex> enc(builder.enc_);
            const std::pair<tax_t, std::string> gs[] = {{2, g2}, {3, g3}, {4, g4a + island + g4b}};
            for(const auto &g: gs)
                enc.for_each([&](u64 km) {
                    const khint_t ki = kh_put(c, db, km, &khr);
                    kh_val(db, ki) = (khr == 0 && kh_val(db, ki) != g.first) ? 1: g.first;
                }, g.second.data(), g.second.size());
        }
        ClassifierGeneric<score::Lex> cl(db, spvec_t{}, k, k, 1, true, false, true, true);
        for(const double d: {0.02, 0.05, 0.10, 0.15})
            for(const unsigned L: {2000u, 5000u, 10000u}) {
                if(isz >= L) continue;
                std::vector<std::string> reads;
                for(unsigned r = 0; r < per_cell; ++r) {
                    std::string rel = mutate(g2.substr(rng() % (g2.size() - L), L - isz), d, rng);
                    const size_t at = rng() % (rel.size() + 1);
                    reads.push_back(rel.substr(0, at) + island + rel.substr(at));
                }
                write_fq(fq, reads);
                for(const Mode &m: MODES) {
                    set_mode(cl, m);
                    const auto calls = classify(cl, taxmap, fq);
                    unsigned tru = 0, lin = 0, wrong = 0, up = 0;
                    for(const tax_t t: calls) {
                        tru += t == 2;
                        lin += t == 2 || t == 10;
                        wrong += t == 4 || t == 20;
                        up += t == 1 || t == 0;
                    }
                    std::printf("%u,%.2f,%u,%s,%zu,%u,%u,%u,%u\n", isz, d, L, m.name, calls.size(), tru, lin, wrong, up);
                }
            }
        if(isz == 1000) {
            // Throughput: error-free reads from all three genomes.
            for(const unsigned L: {150u, 5000u}) {
                const unsigned nr = L == 150 ? 40000: 1200;
                std::vector<std::string> reads;
                const std::string *gs[] = {&g2, &g3, &g4a};
                for(unsigned r = 0; r < nr; ++r) reads.push_back(gs[r % 3]->substr(rng() % (gs[r % 3]->size() - L), L));
                write_fq(fq, reads);
                {   // cost of the KrakenUniq-style counters (-U), pooled
                    ClassifierGeneric<score::Lex> cu(db, spvec_t{}, k, k, 1, true, false, true, true);
                    cu.enable_uniq();
                    classify(cl, taxmap, fq); classify(cu, taxmap, fq);   // warm
                    double tp = 0., tu = 0.;
                    for(int rep = 0; rep < 3; ++rep) {
                        set_mode(cl, MODES[0]);
                        auto t0 = std::chrono::steady_clock::now();
                        classify(cl, taxmap, fq);
                        tp += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
                        t0 = std::chrono::steady_clock::now();
                        classify(cu, taxmap, fq);
                        tu += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
                    }
                    std::fprintf(stderr, "[sweep] -U cost read_bp=%u: %.3f s plain vs %.3f s with counters (+%.1f%%)\n",
                                 L, tp / 3, tu / 3, 100. * (tu - tp) / tp);
                }
                for(const Mode &m: MODES) {
                    set_mode(cl, m);
                    classify(cl, taxmap, fq);   // warm
                    const auto t0 = std::chrono::steady_clock::now();
                    for(int rep = 0; rep < 3; ++rep) classify(cl, taxmap, fq);
                    const double s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() / 3;
                    std::fprintf(stderr, "[sweep] throughput read_bp=%u mode=%s: %.3f s for %u reads (%.2f Mbp/s)\n",
                                 L, m.name, s, nr, nr * double(L) / s / 1e6);
                }
            }
        }
        kh_destroy(c, db);
    }
    kh_destroy(p, taxmap);
}
