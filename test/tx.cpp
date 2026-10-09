#include "test/catch.hpp"

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
        char tmpl[] = "/tmp/bns_tx_XXXXXX";
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
    while(ret.size() < len) ret.push_back(bases[rng() & 3]);
    return ret;
}

std::map<tax_t, tax_t> read_old_to_new(const std::string &path) {
    gzFile fp(gzopen(path.data(), "rb"));
    REQUIRE(fp != nullptr);
    std::map<tax_t, tax_t> ret;
    char buf[256];
    while(gzgets(fp, buf, sizeof(buf))) {
        if(buf[0] == '#') continue;
        tax_t old_id, new_id;
        if(std::sscanf(buf, "%u%u", &old_id, &new_id) == 2) ret[old_id] = new_id;
    }
    gzclose(fp);
    return ret;
}

std::map<std::string, tax_t> read_name_map(const std::string &path) {
    gzFile fp(gzopen(path.data(), "rb"));
    REQUIRE(fp != nullptr);
    std::map<std::string, tax_t> ret;
    char buf[1024];
    while(gzgets(fp, buf, sizeof(buf))) {
        char *tab(std::strchr(buf, '\t'));
        if(!tab) continue;
        *tab = '\0';
        ret[buf] = std::atoi(tab + 1);
    }
    gzclose(fp);
    return ret;
}

std::vector<tax_t> read_new_to_old(const std::string &path) {
    gzFile fp(gzopen(path.data(), "rb"));
    REQUIRE(fp != nullptr);
    std::vector<tax_t> ret;
    tax_t tmp;
    while(gzread(fp, &tmp, sizeof(tmp)) == (int)sizeof(tmp)) ret.push_back(tmp);
    gzclose(fp);
    return ret;
}

} // namespace

TEST_CASE("taxonomy_reformation") {
    TmpDir tmp;
    // NCBI nodes.dmp-style: 1 root; 100, 200 under 1; 101, 102 under 100; 201 under 200.
    const std::string taxpath(tmp.add("tax.txt",
        "100|\t1|\n200|\t1|\n101|\t100|\n102|\t100|\n201|\t200|\n"));
    khash_t(p) *taxmap(build_parent_map(taxpath.data()));
    REQUIRE(taxmap != nullptr);

    // g1 and g2 share taxid 101 (the duplicate-genome case); g3 is 102, g4 is 201.
    const std::string namepath(tmp.add("names.txt", "g1\t101\ng2\t101\ng3\t102\ng4\t201\n"));
    std::mt19937_64 rng(4242);
    std::vector<std::string> paths;
    std::map<std::string, tax_t> expected_old; // keyed by path
    for(const auto &pair: std::vector<std::pair<const char*, tax_t>>{{"g1", 101}, {"g2", 101}, {"g3", 102}, {"g4", 201}}) {
        const std::string path(tmp.add(std::string(pair.first) + ".fa", ">"s + pair.first + "\n" + rand_dna(rng, 200) + "\n"));
        paths.push_back(path);
        expected_old[path] = pair.second;
    }

    TaxonomyReformation tr(namepath.data(), paths, taxmap);
    const std::string o2n(tmp.dir + "/o2n.gz"), n2o(tmp.dir + "/n2o.gz"), nm(tmp.dir + "/nm.gz");
    tr.write_old_to_new(o2n.data());
    tr.write_new_to_old(n2o.data());
    tr.write_name_map(nm.data());
    const auto old_to_new(read_old_to_new(o2n));
    const auto new_to_old(read_new_to_old(n2o));
    const auto name_map(read_name_map(nm));

    SECTION("every genome gets its own node") {
        REQUIRE(name_map.size() == paths.size());
        std::set<tax_t> ids;
        for(const auto &pair: name_map) {
            INFO("genome without a node in name map");
            REQUIRE(expected_old.count(pair.first) == 1);
            ids.insert(pair.second);
        }
        REQUIRE(ids.size() == paths.size()); // distinct nodes
    }
    SECTION("ancestors are present and parents line up") {
        for(const auto old_id: {1u, 100u, 101u, 102u, 200u, 201u}) {
            INFO("expected taxid missing from reformed map");
            REQUIRE(old_to_new.count(old_id) == 1);
        }
        const tax_t n1(old_to_new.at(1)), n100(old_to_new.at(100)), n101(old_to_new.at(101)),
                    n102(old_to_new.at(102)), n200(old_to_new.at(200)), n201(old_to_new.at(201));
        REQUIRE(tr.parent(n100) == n1);
        REQUIRE(tr.parent(n200) == n1);
        REQUIRE(tr.parent(n101) == n100);
        REQUIRE(tr.parent(n102) == n100);
        REQUIRE(tr.parent(n201) == n200);
        REQUIRE(tr.parent(n1) == 0);
    }
    SECTION("duplicate genomes of one taxid are siblings under its node") {
        const tax_t n101(old_to_new.at(101));
        const tax_t g1(name_map.at(paths[0])), g2(name_map.at(paths[1]));
        REQUIRE(g1 != g2);
        REQUIRE(tr.parent(g1) == n101);
        REQUIRE(tr.parent(g2) == n101);
    }
    SECTION("new-to-old round-trips") {
        for(const auto &pair: old_to_new) {
            REQUIRE(pair.second < new_to_old.size());
            REQUIRE(new_to_old[pair.second] == pair.first);
        }
    }
    SECTION("tree is acyclic with a single root") {
        for(const auto &pair: name_map) {
            tax_t cur(pair.second);
            for(unsigned hops(0); cur; ++hops) {
                REQUIRE(hops < 16);
                cur = tr.parent(cur);
            }
        }
    }
    kh_destroy(p, taxmap);
}
