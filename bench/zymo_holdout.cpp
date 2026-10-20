// Read-level adjudication modes on real PacBio HiFi reads of the ZymoBIOMICS
// Gut Microbiome Standard D6331 (SRR13128014, the real-read benchmark of
// Perseus, Nguyen & Schatz 2026): the first 150 MB of the gzipped run
// (~23k reads), against Zymo's D6331 reference genomes (bacteria + archaeon;
// the two yeasts left out to bound memory) plus RefSeq relatives absent from
// the mock (Shigella, Klebsiella, Citrobacter, Enterobacter, Proteus,
// Morganella, Yersinia, Serratia; test/entero/).
// Truth per read (no per-read labels exist): the full-database call under
// hits-only confidence 0.9, kept only if it resolves to a mock species; reads
// without one are excluded.
//  full: every truth read against the full database; "wrong" = any call
//        outside the read's species lineage (e.g. an absent relative).
//  holdouts: a clade's genomes removed; its truth reads classified against the
//        rest. Truth = the deepest ancestor still represented. exact / below
//        (a wrong lower clade, e.g. a wrong family) / above / wrong / unclassified.
// Usage: zymo_holdout <bonsai test dir> > out.csv
#include "classifier.h"

#include <chrono>
#include <cstdio>
#include <cstring>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using namespace bns;

namespace {

struct Genome {const char *file; tax_t taxon; const char *clade;};
#define Z "zymo/refs/D6331.refseq/genomes/"
const Genome GENOMES[] = {
    {Z "Akkermansia_muciniphila.fasta", 239935, ""},
    {Z "Bacteroides_fragilis.fasta", 817, "bacteroidaceae"},
    {Z "Prevotella_corporis.fasta", 28128, ""},
    {Z "Bifidobacterium_adolescentis.fasta", 1680, ""},
    {Z "Clostridioides_difficile.fasta", 1496, ""},
    {Z "Clostridium_perfringens.fasta", 1502, "clostridiaceae"},
    {Z "Faecalibacterium_prausnitzii.fasta", 853, ""},
    {Z "Roseburia_hominis.fasta", 301301, "lachnospiraceae"},
    {Z "Veillonella_rogosae.fasta", 423477, ""},
    {Z "Enterococcus_faecalis.fasta", 1351, "enterococcaceae"},
    {Z "Lactobacillus_fermentum.fasta", 1613, ""},
    {Z "Fusobacterium_nucleatum.fasta", 851, ""},
    {Z "Methanobrevibacter_smithii.fasta", 2173, ""},
    {Z "Salmonella_enterica.fasta", 28901, "entero"},
    {Z "Escherichia_coli_B1109.fasta", 910001, "escherichia"},
    {Z "Escherichia_coli_b2207.fasta", 910002, "escherichia"},
    {Z "Escherichia_coli_B3008.fasta", 910003, "escherichia"},
    {Z "Escherichia_coli_B766.fasta", 910004, "escherichia"},
    {Z "Escherichia_coli_JM109.fasta", 910005, "escherichia"},
    {"entero/shigella.fna.gz", 623, "entero"},
    {"entero/klebsiella.fna.gz", 573, "entero"},
    {"entero/citrobacter.fna.gz", 545, "entero"},
    {"entero/enterobacter.fna.gz", 550, "entero"},
    {"entero/proteus.fna.gz", 584, ""},
    {"entero/morganella.fna.gz", 582, ""},
    {"entero/yersinia.fna.gz", 632, ""},
    {"entero/serratia.fna.gz", 615, ""},
};
#undef Z
// NCBI lineages (main ranks), E. coli strains as local leaves under 562.
const std::pair<tax_t, tax_t> PARENTS[] = {
    {1, 0}, {2, 1}, {2157, 1},
    {74201, 2}, {203494, 74201}, {48461, 203494}, {1647988, 48461}, {239934, 1647988}, {239935, 239934},
    {976, 2}, {200643, 976}, {171549, 200643}, {815, 171549}, {816, 815}, {817, 816},
    {171552, 171549}, {838, 171552}, {28128, 838},
    {201174, 2}, {1760, 201174}, {85004, 1760}, {31953, 85004}, {1678, 31953}, {1680, 1678},
    {1239, 2}, {186801, 1239}, {186802, 186801},
    {186804, 186802}, {1870884, 186804}, {1496, 1870884},
    {31979, 186802}, {1485, 31979}, {1502, 1485},
    {216572, 186802}, {216851, 216572}, {853, 216851},
    {186803, 186802}, {841, 186803}, {301301, 841},
    {909932, 1239}, {1843489, 909932}, {31977, 1843489}, {29465, 31977}, {423477, 29465},
    {91061, 1239}, {186826, 91061}, {81852, 186826}, {1350, 81852}, {1351, 1350},
    {33958, 186826}, {2742598, 33958}, {1613, 2742598},
    {32066, 2}, {203490, 32066}, {203491, 203490}, {203492, 203491}, {848, 203492}, {851, 848},
    {28890, 2157}, {183925, 28890}, {2158, 183925}, {2159, 2158}, {2172, 2159}, {2173, 2172},
    {1224, 2}, {1236, 1224}, {91347, 1236}, {543, 91347},
    {561, 543}, {562, 561}, {910001, 562}, {910002, 562}, {910003, 562}, {910004, 562}, {910005, 562},
    {590, 543}, {28901, 590}, {620, 543}, {623, 620}, {570, 543}, {573, 570},
    {544, 543}, {545, 544}, {547, 543}, {550, 547},
    {1903414, 91347}, {583, 1903414}, {584, 583}, {581, 1903414}, {582, 581},
    {1903411, 91347}, {629, 1903411}, {632, 629}, {613, 1903411}, {615, 613},
};
// Mock species (truth labels); E. coli strains map to 562.
const tax_t MOCK_SPECIES[] = {239935, 817, 28128, 1680, 1496, 1502, 853, 301301, 423477, 1351, 1613, 851, 2173, 28901, 562};

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

std::vector<std::string> read_fasta(const std::string &path) {
    std::vector<std::string> recs;
    std::FILE *p = popen(("gzip -dcf '" + path + "'").data(), "r");
    if(!p) RUNTIME_ERROR("cannot open " + path);
    static char line[1 << 20];
    while(std::fgets(line, sizeof(line), p)) {
        if(line[0] == '>') {recs.emplace_back(); continue;}
        for(char *c = line; *c && *c != '\n' && *c != '\r'; ++c) recs.back().push_back(*c);
    }
    pclose(p);
    return recs;
}
struct Read {std::string name, seq, qual;};
std::vector<Read> read_fastq_slice(const std::string &path) {   // a truncated gzip: drop the partial last record
    std::vector<Read> out;
    std::FILE *p = popen(("gzip -dc '" + path + "' 2>/dev/null").data(), "r");
    static char l[4][1 << 17];
    while(std::fgets(l[0], sizeof(l[0]), p) && std::fgets(l[1], sizeof(l[1]), p)
          && std::fgets(l[2], sizeof(l[2]), p) && std::fgets(l[3], sizeof(l[3]), p)) {
        auto strip = [](char *s) {std::string r(s); while(!r.empty() && (r.back() == '\n' || r.back() == '\r')) r.pop_back(); return r;};
        Read r{strip(l[0]).substr(1), strip(l[1]), strip(l[3])};
        if(r.seq.size() != r.qual.size() || r.seq.empty()) break;
        r.name = r.name.substr(0, r.name.find(' '));
        out.push_back(std::move(r));
    }
    pclose(p);
    return out;
}
bool descends(tax_t a, tax_t root, const khash_t(p) *tm) {   // a strictly below root
    for(unsigned g = 0; a && g < 64; ++g) {
        a = get_parent(tm, a);
        if(a == tax_t(-1)) return false;
        if(a == root) return true;
    }
    return false;
}
tax_t species_of(tax_t t, const khash_t(p) *tm) {
    for(unsigned g = 0; t && t != tax_t(-1) && g < 64; ++g) {
        for(const tax_t s: MOCK_SPECIES) if(s == t) return t;
        t = get_parent(tm, t);
    }
    return 0;
}
template<typename CL>
std::map<std::string, tax_t> classify(const CL &cl, const khash_t(p) *taxmap, const std::string &fq) {
    std::FILE *o = std::tmpfile();
    process_dataset(cl, taxmap, fq.data(), nullptr, o, 1 << 12, 32);
    std::fflush(o); std::rewind(o);
    std::map<std::string, tax_t> ret;
    static char line[1 << 20];
    while(std::fgets(line, sizeof(line), o)) {
        std::istringstream ls(line);
        std::string cu, name;
        tax_t t = 0;
        ls >> cu >> name >> t;
        ret[name] = t;
    }
    std::fclose(o);
    return ret;
}
void write_fq(const std::string &path, const std::vector<const Read *> &reads) {
    std::FILE *a = std::fopen(path.data(), "w");
    for(const Read *r: reads) std::fprintf(a, "@%s\n%s\n+\n%s\n", r->name.data(), r->seq.data(), r->qual.data());
    std::fclose(a);
}
khash_t(c) *build(const std::vector<std::vector<std::string>> &seqs, const std::set<std::string> &removed,
                  const khash_t(p) *taxmap, unsigned k) {
    khash_t(c) *db(kh_init(c));
    kh_resize(c, db, 1u << 27);
    ClassifierGeneric<score::Lex> builder(db, spvec_t{}, k, k, 1, true, false, true, true);
    KmerKeys<score::Lex> enc(builder.enc_);
    int khr;
    for(size_t gi = 0; gi < std::size(GENOMES); ++gi) {
        if(removed.count(GENOMES[gi].clade)) continue;
        const tax_t t = GENOMES[gi].taxon;
        for(const auto &rec: seqs[gi])
            enc.for_each([&](u64 km) {
                const khint_t ki = kh_put(c, db, km, &khr);
                if(khr) kh_val(db, ki) = t;
                else if(kh_val(db, ki) != t) kh_val(db, ki) = lca(taxmap, kh_val(db, ki), t);
            }, rec.data(), rec.size());
    }
    return db;
}

} // namespace

int main(int argc, char **argv) {
    if(argc < 2) {std::fprintf(stderr, "usage: %s <bonsai test dir>\n", *argv); return 1;}
    const std::string dir = argv[1];
    const unsigned k = 31;
    khash_t(p) *taxmap(kh_init(p));
    int khr;
    for(const auto &pc: PARENTS) {const khint_t ki = kh_put(p, taxmap, pc.first, &khr); kh_val(taxmap, ki) = pc.second;}
    std::vector<std::vector<std::string>> seqs;
    for(const auto &g: GENOMES) seqs.push_back(read_fasta(dir + "/" + g.file));
    const auto reads = read_fastq_slice(dir + "/zymo/SRR13128014.head150M.fastq.gz");
    size_t bp = 0;
    for(const auto &r: reads) bp += r.seq.size();
    std::fprintf(stderr, "[zymo] %zu reads, %zu bp; k = %u\n", reads.size(), bp, k);
    char tmpl[] = "/tmp/bns_zymo_XXXXXX";
    const std::string tdir(mkdtemp(tmpl)), fq_all = tdir + "/all.fq", fq_sub = tdir + "/sub.fq";
    std::vector<const Read *> all;
    for(const auto &r: reads) all.push_back(&r);
    write_fq(fq_all, all);
    std::printf("setting,truth_group,mode,reads,exact,below,above,wrong,unclassified\n");

    // Full database: truth labels, then every mode on the truth reads.
    std::map<std::string, tax_t> truth;   // read -> mock species
    {
        auto t0 = std::chrono::steady_clock::now();
        khash_t(c) *db = build(seqs, {}, taxmap, k);
        std::fprintf(stderr, "[zymo] full database: %u k-mers, %.1f s\n", kh_size(db),
                     std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
        ClassifierGeneric<score::Lex> cl(db, spvec_t{}, k, k, 1, true, false, true, true);
        cl.set_confidence(0.9); cl.set_confidence_hits_only(true);
        for(const auto &kv: classify(cl, taxmap, fq_all))
            if(const tax_t s = species_of(kv.second, taxmap); s && (kv.second == s || descends(kv.second, s, taxmap)))
                truth[kv.first] = s;
        std::map<tax_t, size_t> per;
        for(const auto &kv: truth) ++per[kv.second];
        std::fprintf(stderr, "[zymo] truth reads: %zu of %zu;", truth.size(), reads.size());
        for(const auto &kv: per) std::fprintf(stderr, " %u:%zu", kv.first, kv.second);
        std::fprintf(stderr, "\n");
        std::vector<const Read *> sub;
        for(const auto &r: reads) if(truth.count(r.name)) sub.push_back(&r);
        write_fq(fq_sub, sub);
        for(const Mode &m: MODES) {
            cl.set_confidence(m.conf); cl.set_confidence_hits_only(m.hits_only); cl.set_confidence_quality(m.quality);
            const auto t1 = std::chrono::steady_clock::now();
            const auto calls = classify(cl, taxmap, fq_sub);
            const double s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t1).count();
            unsigned ex = 0, below = 0, above = 0, wrong = 0, unc = 0;
            for(const auto &kv: calls) {
                const tax_t T = truth.at(kv.first), C = kv.second;
                if(!C) ++unc;
                else if(C == T) ++ex;
                else if(descends(C, T, taxmap)) ++below;   // a strain of the right species: counted below, not wrong
                else if(descends(T, C, taxmap)) ++above;
                else ++wrong;
            }
            std::printf("full,all,%s,%zu,%u,%u,%u,%u,%u\n", m.name, calls.size(), ex, below, above, wrong, unc);
            std::fprintf(stderr, "[zymo] full %s: %.2f s (%.1f Mbp/s)\n", m.name, s, double(bp) * sub.size() / reads.size() / s / 1e6);
            std::fflush(stdout);
        }
        kh_destroy(c, db);
    }
    // Holdouts: clade removed; truth = deepest ancestor still represented.
    struct Holdout {const char *name; std::set<std::string> removed; std::vector<tax_t> species; tax_t truth;};
    const Holdout HOLDOUTS[] = {
        {"escherichia", {"escherichia"}, {562}, 543},
        {"enterobacteriaceae", {"escherichia", "entero"}, {562, 28901}, 91347},
        {"bacteroidaceae", {"bacteroidaceae"}, {817}, 171549},
        {"lachnospiraceae", {"lachnospiraceae"}, {301301}, 186802},
        {"clostridiaceae", {"clostridiaceae"}, {1502}, 186802},
        {"enterococcaceae", {"enterococcaceae"}, {1351}, 186826},
    };
    for(const Holdout &ho: HOLDOUTS) {
        auto t0 = std::chrono::steady_clock::now();
        khash_t(c) *db = build(seqs, ho.removed, taxmap, k);
        std::vector<const Read *> sub;
        for(const auto &r: reads) {
            auto it = truth.find(r.name);
            if(it == truth.end()) continue;
            for(const tax_t s: ho.species) if(s == it->second) {sub.push_back(&r); break;}
        }
        write_fq(fq_sub, sub);
        std::fprintf(stderr, "[zymo] holdout %s: %u k-mers, %zu reads, built in %.1f s\n", ho.name, kh_size(db), sub.size(),
                     std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
        ClassifierGeneric<score::Lex> cl(db, spvec_t{}, k, k, 1, true, false, true, true);
        for(const Mode &m: MODES) {
            cl.set_confidence(m.conf); cl.set_confidence_hits_only(m.hits_only); cl.set_confidence_quality(m.quality);
            const auto calls = classify(cl, taxmap, fq_sub);
            unsigned ex = 0, below = 0, above = 0, wrong = 0, unc = 0;
            for(const auto &kv: calls) {
                const tax_t C = kv.second;
                if(!C) ++unc;
                else if(C == ho.truth) ++ex;
                else if(descends(C, ho.truth, taxmap)) ++below;
                else if(descends(ho.truth, C, taxmap)) ++above;
                else ++wrong;
            }
            std::printf("%s,%s,%s,%zu,%u,%u,%u,%u,%u\n", ho.name, ho.name, m.name, calls.size(), ex, below, above, wrong, unc);
            std::fflush(stdout);
        }
        kh_destroy(c, db);
    }
    kh_destroy(p, taxmap);
}
