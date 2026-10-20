// Read-level adjudication modes on real genomes with held-out taxa.
// Genomes (bonsai test/ and test/ec/): five E. coli strains, three
// Methanobacterium formicicum strains, Haladaptatus sp. W1, with their NCBI
// lineages (assembly-level leaves get local ids >= 900000 where the strain has
// no taxid of its own). Holdouts:
//   sakai:   E. coli O157:H7 Sakai held out (truth: species 562). A mosaic genome:
//            prophages and pathogenicity islands shared unevenly with other strains.
//   cft073:  E. coli CFT073 held out (truth 562).
//   methano: every Methanobacterium held out (truth: Euryarchaeota 28890, whose only
//            remaining member is Haladaptatus, a different class); real spurious hits.
// Reads are simulated from the held-out genome(s), with a constant Phred
// quality matching the error rate (Q40 error-free): 150 bp error-free, and
// 2/5/10 kb with "hifi" (0.5% subs) or "ont" (3% subs, 1% ins, 1% del) errors.
// For each mode, the call C vs the truth T: exact (C == T), below (C a strict
// descendant of T: a wrong strain or clade), above (strict ancestor), wrong
// (unrelated), unclassified.
// Usage: real_holdout <bonsai test dir> [reads=300] [seed=11] > out.csv
#include "classifier.h"

#include <chrono>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <map>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using namespace bns;

namespace {

struct Genome {const char *file; tax_t taxon; const char *group;};
const Genome GENOMES[] = {
    {"ec/GCF_000005845.2_ASM584v2_genomic.fna.gz", 511145, "ecoli"},   // K-12 MG1655
    {"ec/GCF_000010245.2_ASM1024v1_genomic.fna.gz", 316407, "ecoli"},  // K-12 W3110
    {"ec/GCF_000009565.1_ASM956v1_genomic.fna.gz", 469008, "ecoli"},   // BL21(DE3)
    {"ec/GCF_000007445.1_ASM744v1_genomic.fna.gz", 199310, "cft073"},  // CFT073
    {"ec/GCF_000008865.1_ASM886v1_genomic.fna.gz", 386585, "sakai"},   // O157:H7 Sakai
    {"GCF_000302455.1_ASM30245v1_genomic.fna.gz", 900001, "methano"},  // M. formicicum DSM 3637
    {"GCF_000762265.1_ASM76226v1_genomic.fna.gz", 900002, "methano"},  // M. formicicum BRM9
    {"GCF_000953115.1_DSM1535_genomic.fna.gz", 900003, "methano"},     // M. formicicum DSM 1535
    {"GCF_001723155.1_ASM172315v1_genomic.fna.gz", 900004, "halo"},    // Haladaptatus sp. W1
};
const std::pair<tax_t, tax_t> PARENTS[] = {
    {1, 0}, {2, 1}, {1224, 2}, {1236, 1224}, {91347, 1236}, {543, 91347}, {561, 543}, {562, 561},
    {83333, 562}, {511145, 83333}, {316407, 83333}, {469008, 562}, {199310, 562}, {386585, 562},
    {2157, 1}, {28890, 2157}, {183925, 28890}, {2158, 183925}, {2159, 2158}, {2160, 2159}, {2162, 2160},
    {900001, 2162}, {900002, 2162}, {900003, 2162},
    {183963, 28890}, {2235, 183963}, {367188, 2235}, {900004, 367188},
};
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

std::vector<std::string> read_fasta_gz(const std::string &path) {
    std::vector<std::string> recs;
    std::FILE *p = popen(("gzip -dc '" + path + "'").data(), "r");
    if(!p) RUNTIME_ERROR("cannot open " + path);
    char line[1 << 16];
    while(std::fgets(line, sizeof(line), p)) {
        if(line[0] == '>') {recs.emplace_back(); continue;}
        for(char *c = line; *c && *c != '\n'; ++c) recs.back().push_back(*c);
    }
    pclose(p);
    return recs;
}
std::string mutate(const std::string &s, double sub, double ins, double del, std::mt19937_64 &rng) {
    static const char b[] = "ACGT";
    std::uniform_real_distribution<double> u;
    std::string o;
    o.reserve(s.size() + s.size() / 50);
    for(char c: s) {
        const double x = u(rng);
        if(x < del) continue;
        if(x < del + ins) o.push_back(b[rng() & 3]);
        if(u(rng) < sub) {char n; do n = b[rng() & 3]; while(n == c); c = n;}
        o.push_back(c);
    }
    return o;
}
bool descends(tax_t a, tax_t root, const khash_t(p) *tm) {   // a strictly below root
    for(unsigned g = 0; a && g < 64; ++g) {
        a = get_parent(tm, a);
        if(a == tax_t(-1)) return false;
        if(a == root) return true;
    }
    return false;
}

} // namespace

int main(int argc, char **argv) {
    if(argc < 2) {std::fprintf(stderr, "usage: %s <bonsai test dir> [reads=300] [seed=11]\n", *argv); return 1;}
    const std::string dir = argv[1];
    const unsigned nreads = argc > 2 ? std::atoi(argv[2]): 300;
    const uint64_t seed = argc > 3 ? std::strtoull(argv[3], nullptr, 10): 11;
    const unsigned k = 31;
    std::fprintf(stderr, "[real] seed %llu, %u reads per condition, k = %u\n", (unsigned long long)seed, nreads, k);
    std::mt19937_64 rng(seed);
    khash_t(p) *taxmap(kh_init(p));
    int khr;
    for(const auto &pc: PARENTS) {const khint_t ki = kh_put(p, taxmap, pc.first, &khr); kh_val(taxmap, ki) = pc.second;}
    std::vector<std::vector<std::string>> seqs;
    for(const auto &g: GENOMES) seqs.push_back(read_fasta_gz(dir + "/" + g.file));
    char tmpl[] = "/tmp/bns_real_XXXXXX";
    const std::string fq = std::string(mkdtemp(tmpl)) + "/r.fq";
    std::printf("holdout,read_bp,errors,mode,reads,exact,below,above,wrong,unclassified\n");
    struct Holdout {const char *group; tax_t truth;};
    for(const Holdout ho: {Holdout{"sakai", 562}, Holdout{"cft073", 562}, Holdout{"methano", 28890}}) {
        const auto t0 = std::chrono::steady_clock::now();
        khash_t(c) *db(kh_init(c));
        size_t ninserted = 0;
        {
            ClassifierGeneric<score::Lex> builder(db, spvec_t{}, k, k, 1, true, false, true, true);
            KmerKeys<score::Lex> enc(builder.enc_);
            for(size_t gi = 0; gi < std::size(GENOMES); ++gi) {
                if(!std::strcmp(GENOMES[gi].group, ho.group)) continue;
                const tax_t t = GENOMES[gi].taxon;
                for(const auto &rec: seqs[gi])
                    enc.for_each([&](u64 km) {
                        const khint_t ki = kh_put(c, db, km, &khr);
                        if(khr) kh_val(db, ki) = t, ++ninserted;
                        else if(kh_val(db, ki) != t) kh_val(db, ki) = lca(taxmap, kh_val(db, ki), t);
                    }, rec.data(), rec.size());
            }
        }
        std::fprintf(stderr, "[real] holdout %s: %zu k-mers, built in %.1f s\n", ho.group, ninserted,
                     std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
        // Held-out sequence, records weighted by length
        std::vector<const std::string *> src;
        for(size_t gi = 0; gi < std::size(GENOMES); ++gi)
            if(!std::strcmp(GENOMES[gi].group, ho.group)) for(const auto &r: seqs[gi]) src.push_back(&r);
        size_t total = 0;
        for(auto *r: src) total += r->size();
        ClassifierGeneric<score::Lex> cl(db, spvec_t{}, k, k, 1, true, false, true, true);
        struct Cond {unsigned L; const char *name; double sub, ins, del;};
        for(const Cond cd: {Cond{150, "none", 0, 0, 0}, Cond{2000, "hifi", 0.005, 0, 0}, Cond{5000, "hifi", 0.005, 0, 0},
                            Cond{10000, "hifi", 0.005, 0, 0}, Cond{2000, "ont", 0.03, 0.01, 0.01},
                            Cond{5000, "ont", 0.03, 0.01, 0.01}, Cond{10000, "ont", 0.03, 0.01, 0.01}}) {
            {
                std::FILE *a = std::fopen(fq.data(), "w");
                for(unsigned r = 0; r < nreads;) {
                    size_t off = rng() % total, ri = 0;
                    while(off >= src[ri]->size()) off -= src[ri++]->size();
                    if(src[ri]->size() < cd.L + 1 || off + cd.L > src[ri]->size()) continue;
                    const std::string s = mutate(src[ri]->substr(off, cd.L), cd.sub, cd.ins, cd.del, rng);
                    // constant Phred quality matching the simulated per-base error rate (Q40 when error-free)
                    const double perr = cd.sub + cd.ins + cd.del;
                    const int q = perr > 0 ? std::min(40, int(std::lround(-10 * std::log10(perr)))): 40;
                    std::fprintf(a, "@r%u\n%s\n+\n%s\n", r++, s.data(), std::string(s.size(), char(33 + q)).data());
                }
                std::fclose(a);
            }
            for(const Mode &m: MODES) {
                cl.set_confidence(m.conf); cl.set_confidence_hits_only(m.hits_only); cl.set_confidence_quality(m.quality);
                std::FILE *o = std::tmpfile();
                process_dataset(cl, taxmap, fq.data(), nullptr, o, 1 << 14, 32);
                std::fflush(o); std::rewind(o);
                unsigned n = 0, ex = 0, below = 0, above = 0, wrong = 0, unc = 0;
                char line[1 << 16];
                while(std::fgets(line, sizeof(line), o)) {
                    std::istringstream ls(line);
                    std::string cu, name;
                    tax_t t = 0;
                    ls >> cu >> name >> t;
                    ++n;
                    if(!t) ++unc;
                    else if(t == ho.truth) ++ex;
                    else if(descends(t, ho.truth, taxmap)) ++below;
                    else if(descends(ho.truth, t, taxmap)) ++above;
                    else ++wrong;
                }
                std::fclose(o);
                std::printf("%s,%u,%s,%s,%u,%u,%u,%u,%u,%u\n", ho.group, cd.L, cd.name, m.name, n, ex, below, above, wrong, unc);
                std::fflush(stdout);
            }
        }
        kh_destroy(c, db);
    }
    kh_destroy(p, taxmap);
}
