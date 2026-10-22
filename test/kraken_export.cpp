#include "test/catch.hpp"

#include "kraken_export.h"
#include "tx.h"
#include "util.h"

#include <map>
#include <random>
#include <string>
#include <vector>

using namespace bns;

namespace {

struct TmpDir {
    std::string dir;
    std::vector<std::string> files;
    TmpDir(): dir([] {
        char tmpl[] = "/tmp/bns_kraken_XXXXXX";
        const char *ret(mkdtemp(tmpl));
        REQUIRE(ret != nullptr);
        return std::string(ret);
    }()) {}
    ~TmpDir() {
        for(const auto &f: files) std::remove(f.data());
        std::remove(dir.data());
    }
    std::string add(const std::string &name, const std::string &content) {
        files.push_back(dir + "/" + name);
        std::FILE *fp(std::fopen(files.back().data(), "w"));
        REQUIRE(fp != nullptr);
        std::fputs(content.data(), fp);
        std::fclose(fp);
        return files.back();
    }
};

std::string rand_dna(std::mt19937_64 &rng, size_t len) {
    static const char bases[] = "ACGT";
    std::string ret;
    ret.reserve(len);
    while(ret.size() < len) ret.push_back(bases[(rng() >> 33) & 3]);
    return ret;
}

// Parse nodes.dmp the way kraken2's taxonomy.cc does: strip the trailing
// two characters, split on "\t|\t", take id, parent, rank.
std::map<tax_t, std::pair<tax_t, std::string>> parse_nodes_dmp(const std::string &path) {
    std::FILE *fp(std::fopen(path.data(), "r"));
    REQUIRE(fp != nullptr);
    std::map<tax_t, std::pair<tax_t, std::string>> ret;
    char *line(nullptr);
    size_t cap(0);
    ssize_t len;
    while((len = getline(&line, &cap, fp)) >= 0) {
        std::string l(line);
        while(!l.empty() && (l.back() == '\n' || l.back() == '\r')) l.pop_back();
        REQUIRE(l.size() >= 2);
        l.resize(l.size() - 2); // discard trailing "\t|"
        std::vector<std::string> fields;
        size_t pos(0);
        while(true) {
            const auto npos(l.find("\t|\t", pos));
            if(npos == std::string::npos) {fields.push_back(l.substr(pos)); break;}
            fields.push_back(l.substr(pos, npos - pos));
            pos = npos + 3;
        }
        REQUIRE(fields.size() >= 3);
        ret[std::atoi(fields[0].data())] = {std::atoi(fields[1].data()), fields[2]};
    }
    std::free(line);
    std::fclose(fp);
    return ret;
}

std::map<tax_t, std::string> parse_names_dmp(const std::string &path) {
    std::FILE *fp(std::fopen(path.data(), "r"));
    REQUIRE(fp != nullptr);
    std::map<tax_t, std::string> ret;
    char buf[4096];
    while(std::fgets(buf, sizeof(buf), fp)) {
        std::string l(buf);
        while(!l.empty() && (l.back() == '\n' || l.back() == '\r')) l.pop_back();
        l.resize(l.size() - 2);
        tax_t id(std::atoi(l.data()));
        const auto p1(l.find("\t|\t"));
        REQUIRE(p1 != std::string::npos);
        const auto name(l.substr(p1 + 3, l.find("\t|\t", p1 + 3) - p1 - 3));
        REQUIRE(l.find("scientific name") != std::string::npos);
        ret[id] = name;
    }
    std::fclose(fp);
    return ret;
}

std::vector<tax_t> ancestors(const std::map<tax_t, std::pair<tax_t, std::string>> &nodes, tax_t t) {
    std::vector<tax_t> ret;
    while(true) {
        ret.push_back(t);
        const auto it(nodes.find(t));
        REQUIRE(it != nodes.end());
        if(it->second.first == t) break; // root is its own parent
        t = it->second.first;
        REQUIRE(ret.size() < 64);
    }
    return ret;
}

} // namespace

TEST_CASE("kraken_export") {
    TmpDir tmp;
    std::mt19937_64 rng(555);
    enum {NGEN = 8};
    const std::string taxpath(tmp.add("tax.txt", "100|\t1|\n101|\t1|\n102|\t1|\n103|\t1|\n"));
    khash_t(p) *taxmap(build_parent_map(taxpath.data()));
    std::string names;
    std::vector<std::string> paths, seqids;
    const std::string species_blocks[4] = {rand_dna(rng, 50), rand_dna(rng, 50), rand_dna(rng, 50), rand_dna(rng, 50)};
    const std::string module_block(rand_dna(rng, 50));
    for(unsigned g(0); g < NGEN; ++g) {
        const unsigned sp(g >> 1);
        const std::string name("g" + std::to_string(g));
        std::string seq(rand_dna(rng, 400) + species_blocks[sp]);
        if(sp == 0 || sp == 2) seq += module_block;
        paths.push_back(tmp.add(name + ".fa", ">" + name + "\n" + seq + "\n"));
        seqids.push_back(name);
        names += name + "\t" + std::to_string(100 + sp) + "\n";
    }
    const std::string namepath(tmp.add("names.txt", names));
    TaxonomyReformation tr(namepath.data(), paths, taxmap);
    const std::string nm(tmp.dir + "/nm.gz");
    tr.write_name_map(nm.data());
    gzFile fp(gzopen(nm.data(), "rb"));
    REQUIRE(fp != nullptr);
    std::map<std::string, tax_t> path_to_new;
    char buf[2048];
    while(gzgets(fp, buf, sizeof(buf))) {
        char *tab(std::strchr(buf, '\t'));
        REQUIRE(tab != nullptr);
        *tab = '\0';
        path_to_new[buf] = std::atoi(tab + 1);
    }
    gzclose(fp);
    std::vector<tax_t> new_gt;
    for(const auto &path: paths) new_gt.push_back(path_to_new.at(path));
    PlacementTree tree(tr.pmap(), new_gt);
    Spacer spacer(17, 17);
    KmerClasses kc;
    kc.fill<score::Lex>(paths, spacer);
    const GreedyResult gr(greedy_synthetic_nodes(tree, kc, 4));
    REQUIRE(gr.added.size() == 1); // the module node
    REQUIRE(tree.total_delta(kc) == 0);

    const std::string outdir(tmp.dir + "/export");
    ::mkdir(outdir.data(), 0755);
    const auto ret(KrakenExporter::write(tree, seqids, paths, outdir.data(), "grafted_"));
    REQUIRE(ret.seqid2taxid.size() == NGEN);
    REQUIRE(ret.ngrafted == 1);

    const auto nodes(parse_nodes_dmp(outdir + "/nodes.dmp"));
    const auto names_out(parse_names_dmp(outdir + "/names.dmp"));
    SECTION("nodes.dmp structure") {
        REQUIRE(nodes.size() == names_out.size());
        // Root is its own parent and keeps taxonomy id 1.
        const auto root_it(nodes.find(1));
        REQUIRE(root_it != nodes.end());
        REQUIRE(root_it->second.first == 1);
        // Reformed taxonomy ids are preserved in the export.
        for(const auto id: new_gt) REQUIRE(nodes.count(id) == 1);
        // Every node reaches the root.
        for(const auto &pair: nodes) {
            const auto anc(ancestors(nodes, pair.first));
            REQUIRE(anc.back() == 1);
        }
    }
    SECTION("grafted node is marked and named") {
        tax_t grafted(0);
        for(const auto &pair: nodes)
            if(names_out.at(pair.first).rfind("grafted_", 0) == 0) grafted = pair.first;
        REQUIRE(grafted != 0);
        REQUIRE(nodes.at(grafted).second == "no rank");
        REQUIRE(nodes.count(nodes.at(grafted).first) == 1); // parent exists
        // It must sit above the leaves of the four module genomes but not
        // the others: the export preserves the tree structure.
        std::map<std::string, tax_t> leaf_of;
        for(const auto &pair: ret.seqid2taxid) leaf_of[pair.first] = pair.second;
        for(unsigned g(0); g < NGEN; ++g) {
            const auto anc(ancestors(nodes, leaf_of["g" + std::to_string(g)]));
            REQUIRE(std::find(anc.begin(), anc.end(), new_gt[g]) != anc.end()); // reformed per-genome node on the path
            const bool under(std::find(anc.begin(), anc.end(), grafted) != anc.end());
            REQUIRE(under == (g <= 1 || g == 4 || g == 5)); // species 0 and 2 hold the module
        }
    }
    SECTION("seqid2taxid.map and libraries") {
        std::FILE *fp(std::fopen((outdir + "/seqid2taxid.map").data(), "r"));
        REQUIRE(fp != nullptr);
        char line[4096];
        unsigned n(0);
        while(std::fgets(line, sizeof(line), fp)) {
            const auto it(ret.seqid2taxid[n]);
            std::string expect(it.first + "\t" + std::to_string(it.second) + "\n");
            REQUIRE(line == expect);
            ++n;
        }
        std::fclose(fp);
        REQUIRE(n == NGEN);
        // Kraken 2's build_db keys its map on the full rewritten header
        // token, not the bare seqid.
        std::FILE *f3(std::fopen((outdir + "/seqid2taxid_kraken2.map").data(), "r"));
        REQUIRE(f3 != nullptr);
        REQUIRE(std::fgets(line, sizeof(line), f3) != nullptr);
        {
            const auto it(ret.seqid2taxid[0]);
            std::string expect(it.first + "|kraken:taxid|" + std::to_string(it.second) + "\t" + std::to_string(it.second) + "\n");
            REQUIRE(line == expect);
        }
        std::fclose(f3);
        // Kraken 2 library headers carry the explicit taxid.
        std::FILE *f2(std::fopen((outdir + "/library_kraken2/g0.fna").data(), "r"));
        REQUIRE(f2 != nullptr);
        char h[4096];
        REQUIRE(std::fgets(h, sizeof(h), f2) != nullptr);
        std::string got(h);
        got.pop_back(); // newline
        REQUIRE(got == (">" + seqids[0] + "|kraken:taxid|" + std::to_string(ret.seqid2taxid[0].second)));
        std::fclose(f2);
        // KrakenUniq library is a verbatim copy.
        std::FILE *a(std::fopen(paths[0].data(), "r")), *b(std::fopen((outdir + "/library/g0.fna").data(), "r"));
        REQUIRE(a != nullptr);
        REQUIRE(b != nullptr);
        int ca(0), cb(0);
        bool same(true);
        for(;;) {
            ca = std::fgetc(a);
            cb = std::fgetc(b);
            if(ca != cb) {same = false; break;}
            if(ca == EOF) break;
        }
        std::fclose(a);
        std::fclose(b);
        REQUIRE(same);
    }
    kh_destroy(p, taxmap);
}

TEST_CASE("kraken_export dedupes colliding library stems") {
    TmpDir tmp;
    std::mt19937_64 rng(777);
    ::mkdir((tmp.dir + "/a").data(), 0755);
    ::mkdir((tmp.dir + "/b").data(), 0755);
    const std::string taxpath(tmp.add("tax.txt", "100|\t1|\n101|\t1|\n"));
    khash_t(p) *taxmap(build_parent_map(taxpath.data()));
    std::vector<std::string> paths, seqids;
    std::vector<tax_t> genome_taxa;
    const char *subs[] {"/a", "/b"};
    const tax_t tids[] {100, 101};
    for(unsigned i(0); i < 2; ++i) {
        const std::string name(std::string(subs[i] + 1) + "genome"); // agenome, bgenome
        paths.push_back(tmp.add(std::string(subs[i] + 1) + "/same.fa", ">" + name + "\n" + rand_dna(rng, 200) + "\n"));
        seqids.push_back(name);
        genome_taxa.push_back(tids[i]);
    }
    PlacementTree tree(taxmap, genome_taxa);
    const std::string outdir(tmp.dir + "/export2");
    ::mkdir(outdir.data(), 0755);
    const auto ret(KrakenExporter::write(tree, seqids, paths, outdir.data()));
    REQUIRE(ret.seqid2taxid.size() == 2);
    // Both library copies must exist under distinct stems; one genome's
    // fasta must not have overwritten the other's.
    std::FILE *fp(std::fopen((outdir + "/library/same.fna").data(), "r"));
    REQUIRE(fp != nullptr);
    std::fclose(fp);
    fp = std::fopen((outdir + "/library/same_2.fna").data(), "r");
    REQUIRE(fp != nullptr);
    std::fclose(fp);
    fp = std::fopen((outdir + "/library_kraken2/same.fna").data(), "r");
    REQUIRE(fp != nullptr);
    std::fclose(fp);
    fp = std::fopen((outdir + "/library_kraken2/same_2.fna").data(), "r");
    REQUIRE(fp != nullptr);
    std::fclose(fp);
    kh_destroy(p, taxmap);
}
