// Six-frame translation (sixframe.h) and translated protein databases
// (KeySpec{PROTEIN20, translate}):
//  - all 64 codons vs an independently written NCBI table 1; frames +1..+3 and
//    -1..-3 on a known sequence; stops split peptides; N-codons give X;
//  - a compact protein database built from 6-frame-translated genomes, and
//    one built from protein FASTA references (frame +1 of each genome),
//    both classify error-free nucleotide reads to their own genome's taxon;
//    the KeySpec round-trips through the database file;
//  - reads with SYNONYMOUS third-position changes in the genome's frame keep
//    their protein k-mers and still classify through the protein database,
//    while the DNA (k = 31) database loses most of them. Two genome models:
//    an ideal open reading frame (stop-free frame +1: every synonymous read
//    classifies) and plain random DNA with stops wherever they fall (a read
//    is lost only if its frame has no stop-free run of k residues).
#include "test/catch.hpp"

#include "classifier.h"
#include "feature_min.h"
#include "kmermap.h"
#include "sixframe.h"

#include <cstdio>
#include <map>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using namespace bns;

namespace {

// random DNA; coding = true draws codons rejecting stops in frame +1 (an ideal
// open reading frame), coding = false leaves stops wherever they fall (as
// introns, pseudogenes, frameshifts and sequencing errors produce them)
std::string random_genome(size_t ncodons, bool coding, std::mt19937_64 &rng) {
    static const char b[] = "ACGT";
    std::string s;
    s.reserve(3 * ncodons);
    while(s.size() < 3 * ncodons) {
        char c[4] = {b[rng() & 3], b[rng() & 3], b[rng() & 3], 0};
        if(coding && sixframe::translate_codon(c) == '*') continue;
        s.append(c, 3);
    }
    return s;
}
void put(const std::string &path, const std::string &content) {
    std::FILE *fp = std::fopen(path.data(), "w");
    REQUIRE(fp != nullptr);
    std::fputs(content.data(), fp);
    std::fclose(fp);
}
template<typename CL>
std::map<std::string, tax_t> classify(const CL &cl, const khash_t(p) *taxmap, const std::string &fq) {
    std::FILE *o = std::tmpfile();
    process_dataset(cl, taxmap, fq.data(), nullptr, o, 256, 32);
    std::fflush(o); std::rewind(o);
    std::string out;
    char buf[1 << 14];
    size_t n;
    while((n = std::fread(buf, 1, sizeof(buf), o)) > 0) out.append(buf, n);
    std::fclose(o);
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

TEST_CASE("six-frame translation: codon table, frames, stops, ambiguous codons") {
    // NCBI translation table 1, written codon by codon
    const std::map<std::string, char> table1 = {
        {"TTT",'F'},{"TTC",'F'},{"TTA",'L'},{"TTG",'L'},{"CTT",'L'},{"CTC",'L'},{"CTA",'L'},{"CTG",'L'},
        {"ATT",'I'},{"ATC",'I'},{"ATA",'I'},{"ATG",'M'},{"GTT",'V'},{"GTC",'V'},{"GTA",'V'},{"GTG",'V'},
        {"TCT",'S'},{"TCC",'S'},{"TCA",'S'},{"TCG",'S'},{"CCT",'P'},{"CCC",'P'},{"CCA",'P'},{"CCG",'P'},
        {"ACT",'T'},{"ACC",'T'},{"ACA",'T'},{"ACG",'T'},{"GCT",'A'},{"GCC",'A'},{"GCA",'A'},{"GCG",'A'},
        {"TAT",'Y'},{"TAC",'Y'},{"TAA",'*'},{"TAG",'*'},{"CAT",'H'},{"CAC",'H'},{"CAA",'Q'},{"CAG",'Q'},
        {"AAT",'N'},{"AAC",'N'},{"AAA",'K'},{"AAG",'K'},{"GAT",'D'},{"GAC",'D'},{"GAA",'E'},{"GAG",'E'},
        {"TGT",'C'},{"TGC",'C'},{"TGA",'*'},{"TGG",'W'},{"CGT",'R'},{"CGC",'R'},{"CGA",'R'},{"CGG",'R'},
        {"AGT",'S'},{"AGC",'S'},{"AGA",'R'},{"AGG",'R'},{"GGT",'G'},{"GGC",'G'},{"GGA",'G'},{"GGG",'G'}};
    REQUIRE(table1.size() == 64);
    for(const auto &kv: table1) {
        REQUIRE(sixframe::translate_codon(kv.first.data()) == kv.second);
        std::string lower(kv.first), rna(kv.first);
        for(auto &c: lower) c = char(std::tolower(c));
        for(auto &c: rna) if(c == 'T') c = 'U';
        REQUIRE(sixframe::translate_codon(lower.data()) == kv.second);
        REQUIRE(sixframe::translate_codon(rna.data()) == kv.second);
    }
    REQUIRE(sixframe::translate_codon("ANG") == 'X');
    // ATG GCC TAA AAA CCC: +1 = M A * K P -> peptides "MA", "KP"
    const std::string s = "ATGGCCTAAAAACCC";
    std::vector<std::pair<int, std::string>> peps;
    sixframe::for_each_peptide(s.data(), s.size(), [&](const char *p, size_t n, int frame) {peps.emplace_back(frame, std::string(p, n));});
    std::map<int, std::vector<std::string>> by;
    for(const auto &pf: peps) by[pf.first].push_back(pf.second);
    REQUIRE(by[1] == std::vector<std::string>{"MA", "KP"});
    // +2 from offset 1: TGG CCT AAA AAC = W P K N
    REQUIRE(by[2] == std::vector<std::string>{"WPKN"});
    // reverse complement GGGTTTTTAGGCCAT: -1 = GGG TTT TTA GGC CAT = G F L G H
    REQUIRE(by[-1] == std::vector<std::string>{"GFLGH"});
    REQUIRE(by.count(3));
    REQUIRE(by.count(-2));
    REQUIRE(by.count(-3));
}

static void run_protein_case(bool coding) {
    char tmpl[] = "/tmp/bns_6f_XXXXXX";
    const std::string dir(mkdtemp(tmpl));
    std::mt19937_64 rng(31);
    std::vector<std::string> genomes, dna_fns, prot_fns;
    std::string seq2tax;
    for(int g = 0; g < 4; ++g) {
        genomes.push_back(random_genome(10000, coding, rng));
        const std::string name = "g" + std::to_string(g);
        put(dir + "/" + name + ".fa", ">" + name + "\n" + genomes.back() + "\n");
        dna_fns.push_back(dir + "/" + name + ".fa");
        std::string prot;
        sixframe::translate_frame(genomes.back().data(), genomes.back().size(), 0, prot);   // frame +1 as "protein references"
        put(dir + "/" + name + ".faa", ">" + name + "\n" + prot + "\n");
        prot_fns.push_back(dir + "/" + name + ".faa");
        seq2tax += name + "\t" + std::to_string(2 + g) + "\n";
    }
    put(dir + "/seq2tax.txt", seq2tax);
    put(dir + "/tax.txt", "2\t|\t1\n3\t|\t1\n4\t|\t1\n5\t|\t1\n");
    khash_t(p) *taxmap(build_parent_map((dir + "/tax.txt").data()));
    // reads: frame-aligned 150 bp substrings (start = 3j in genome coordinates), clean
    // and with synonymous third-position changes at four-fold degenerate codons
    const std::string fq = dir + "/r.fq";
    std::map<std::string, tax_t> truth;
    {
        std::FILE *a = std::fopen(fq.data(), "w");
        for(int r = 0; r < 400; ++r) {
            const int g = r % 4;
            const size_t start = 3 * (rng() % 9900);
            std::string s = genomes[g].substr(start, 150);
            const bool syn = r >= 200;
            if(syn) {
                for(size_t c = 0; c + 3 <= s.size(); c += 3) {
                    const std::string two = s.substr(c, 2);
                    const bool fourfold = two == "CT" || two == "GT" || two == "TC" || two == "CC" || two == "AC" || two == "GC" || two == "CG" || two == "GG";
                    if(fourfold) {char nb; do nb = "ACGT"[rng() & 3]; while(nb == s[c + 2]); s[c + 2] = nb;}
                }
            }
            const std::string name = (syn ? "syn": "clean") + std::to_string(r);
            truth[name] = tax_t(2 + g);
            std::fprintf(a, "@%s\n%s\n+\n%s\n", name.data(), s.data(), std::string(s.size(), 'I').data());
        }
        std::fclose(a);
    }
    const unsigned k = 12;                     // residues
    Spacer sp(k, k);
    const KeySpec tspec{PROTEIN20, true, 0};
    const unsigned kb = KmerKeys<>::key_bits_for(k, PROTEIN20);
    REQUIRE(kb == 52);                         // bit width of 20^12 - 1
    auto tdb = compact_lca_map<score::Lex>(dna_fns, taxmap, (dir + "/seq2tax.txt").data(), sp, 2, false, 1 << 16, kb, tspec);
    auto pdb = compact_lca_map<score::Lex>(prot_fns, taxmap, (dir + "/seq2tax.txt").data(), sp, 2, false, 1 << 16, kb, KeySpec{PROTEIN20, false, 0});
    REQUIRE(tdb->size() > pdb->size());        // six frames vs one
    // file round trip keeps the key encoding; classify translates reads
    const std::string cdb = dir + "/prot.cdb";
    CompactDatabase::write(cdb.data(), k, k, spvec_t(k - 1, 0), *tdb, tspec);
    CompactDatabase cd(cdb.data());
    REQUIRE(cd.rht_ == PROTEIN20);
    REQUIRE(cd.translate_);
    REQUIRE(cd.keyspec().translate);
    for(const auto *m: {cd.db_.get(), pdb.get()}) {
        ClassifierGeneric<score::Lex, CompactKmerMap> cl(m, spvec_t{}, k, k, 2, true, false, true, false, cd.keyspec());
        REQUIRE(cl.enc_.translated());
        const auto got = classify(cl, taxmap, fq);
        size_t clean_ok = 0, syn_ok = 0;
        for(const auto &kv: truth) {
            const bool ok = got.at(kv.first) == kv.second;
            if(kv.first.rfind("clean", 0) == 0) clean_ok += ok; else syn_ok += ok;
        }
        REQUIRE(clean_ok == 200);
        if(coding) REQUIRE(syn_ok == 200);     // the genome's frame keeps every residue
        else       REQUIRE(syn_ok >= 190);     // stops in the frame can leave < k residues between them
        std::fprintf(stderr, "[%s genomes] protein DB (%s refs, %zu keys): clean %zu/200, synonymous-mutated %zu/200 to own taxon\n",
                     coding ? "coding": "random", m == pdb.get() ? "protein FASTA": "6-frame translated", m->size(), clean_ok, syn_ok);
    }
    // the DNA database loses most synonymous-mutated reads
    {
        Spacer dsp(31, 31);
        auto ddb = compact_lca_map<score::Lex>(dna_fns, taxmap, (dir + "/seq2tax.txt").data(), dsp, 2, true, 1 << 16, 62);
        ClassifierGeneric<score::Lex, CompactKmerMap> dl(ddb.get(), spvec_t{}, 31, 31, 2, true, false, true, true);
        const auto got = classify(dl, taxmap, fq);
        size_t syn_ok = 0, clean_ok = 0;
        for(const auto &kv: truth) {
            const bool ok = got.at(kv.first) == kv.second;
            if(kv.first.rfind("clean", 0) == 0) clean_ok += ok; else syn_ok += ok;
        }
        REQUIRE(clean_ok == 200);
        REQUIRE(syn_ok < 100);
        std::fprintf(stderr, "[%s genomes] DNA DB (k=31): clean %zu/200, synonymous-mutated %zu/200 to own taxon\n",
                     coding ? "coding": "random", clean_ok, syn_ok);
    }
    kh_destroy(p, taxmap);
}

TEST_CASE("translated protein databases classify nucleotide reads; synonymous changes survive") {
    run_protein_case(true);    // ideal open reading frames
    run_protein_case(false);   // stops wherever they fall
}
