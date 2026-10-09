#include "test/catch.hpp"

#include "flextree.h"
#include "util.h"

#include <algorithm>
#include <functional>
#include <random>
#include <string>
#include <vector>
#include <sys/stat.h>

using namespace bns;

namespace {

// Deterministic random DNA string.
std::string rand_dna(std::mt19937_64 &rng, size_t len) {
    static const char bases[] = "ACGT";
    std::string ret;
    ret.reserve(len);
    while(ret.size() < len) ret.push_back(bases[rng() & 3]);
    return ret;
}

std::string write_file(const std::string &path, const std::string &content) {
    std::FILE *fp(std::fopen(path.data(), "w"));
    REQUIRE(fp != nullptr);
    std::fputs(content.data(), fp);
    std::fclose(fp);
    return path;
}

struct TmpGenomes {
    std::string dir;
    std::vector<std::string> paths;
    TmpGenomes(): dir([] {
        char tmpl[] = "/tmp/bns_flextree_XXXXXX";
        const char *ret(mkdtemp(tmpl));
        REQUIRE(ret != nullptr);
        return std::string(ret);
    }()) {}
    ~TmpGenomes() {
        for(const auto &path: paths) std::remove(path.data());
        std::remove(dir.data());
    }
};

/*
Synthetic taxonomy:
  1 (root)
  +- 100
  |   +- 101, 102  (share kmer X)
  |   +- 103, 104  (share kmer Y)
  +- 200
      +- 201, 202  (share kmer Z)
*/
std::string make_taxonomy_file(const std::string &dir) {
    return write_file(dir + "/tax.txt",
        "100|1\n200|1\n101|100\n102|100\n103|100\n104|100\n201|200\n202|200\n");
}

struct EmittedNode {
    tax_t taxid, parent;
    u64 score;
};

std::vector<EmittedNode> parse_emitted(std::FILE *fp) {
    std::rewind(fp);
    std::vector<EmittedNode> ret;
    char *line(nullptr);
    size_t cap(0);
    ssize_t len;
    while((len = getline(&line, &cap, fp)) >= 0) {
        if(line[0] == '#') continue;
        EmittedNode node;
        char *saveptr(nullptr);
        node.taxid = std::strtoul(strtok_r(line, "\t", &saveptr), nullptr, 10);
        node.score = std::strtoull(strtok_r(nullptr, "\t", &saveptr), nullptr, 10);
        node.parent = std::strtoul(strtok_r(nullptr, "\t", &saveptr), nullptr, 10);
        ret.push_back(node);
    }
    std::free(line);
    return ret;
}

} // namespace

TEST_CASE("flextree_greedy_placement") {
    TmpGenomes tmp;
    const std::string taxpath(make_taxonomy_file(tmp.dir));
    khash_t(p) *taxmap(build_parent_map(taxpath.data()));
    REQUIRE(taxmap != nullptr);

    const unsigned k(11);
    Spacer sp(k, k);
    std::mt19937_64 rng(1337);
    const std::string x("ACGTTGCAAGT"), y("TTGACCGATCG"), z("GGCATAGTCAA");
    auto genome = [&](const std::string &special) {
        return rand_dna(rng, 100) + special + rand_dna(rng, 100);
    };
    std::unordered_map<tax_t, std::forward_list<std::string>> tpm;
    auto add_genome = [&](tax_t tax, const std::string &seq) {
        const std::string path(tmp.dir + "/g" + std::to_string(tax) + ".fa");
        write_file(path, ">g" + std::to_string(tax) + "\n" + seq + "\n");
        tmp.paths.push_back(path);
        tpm[tax].push_front(path);
    };
    add_genome(101, genome(x));
    add_genome(102, genome(x));
    add_genome(103, genome(y));
    add_genome(104, genome(y));
    add_genome(201, genome(z));
    add_genome(202, genome(z));

    SECTION("emissions are highest-score-first, recorded, and consistent with output") {
        FMEmitter fme(taxmap, tpm, true, 1 << 8, 11); // 11 - 9 taxa = 2 additions
        std::vector<tax_t> c100{101, 102, 103, 104}, c200{201, 202};
        fme.process_subtree(100, c100.begin(), c100.end(), sp);
        fme.process_subtree(200, c200.begin(), c200.end(), sp);
        std::FILE *fp(tmpfile());
        REQUIRE(fp != nullptr);
        ks::string ks;
        fme.run_collapse(1000, ks, fp);

        const auto &added(fme.added_nodes());
        REQUIRE(added.size() == 2);
        // Only k-mers shared by a proper subset of the subtree survive bitmap_t's
        // filter, so each candidate is a shared pair: score (4 clades - 2 present) * n.
        std::vector<u64> scores;
        for(const auto &pair: added) {
            scores.push_back(pair.second);
            REQUIRE(pair.first == 100);
        }
        REQUIRE(std::is_sorted(scores.begin(), scores.end(), std::greater<u64>()));
        REQUIRE(scores == std::vector<u64>{8, 2}); // (4 - 2) * n: n = 4 shared k-mers (y), then 1 (x)

        const auto emitted(parse_emitted(fp));
        std::fclose(fp);
        REQUIRE(emitted.size() == added.size());
        for(size_t i(0); i < emitted.size(); ++i) {
            REQUIRE(emitted[i].taxid == 1001 + i); // sequential from maxtax + 1
            REQUIRE(emitted[i].parent == added[i].first);
            REQUIRE(emitted[i].score == added[i].second);
        }
    }

    SECTION("requests beyond available candidates stop cleanly once all nodes are added") {
        FMEmitter fme(taxmap, tpm, true, 1 << 8, 100); // 91 additions requested
        std::vector<tax_t> c100{101, 102, 103, 104}, c200{201, 202};
        fme.process_subtree(100, c100.begin(), c100.end(), sp);
        fme.process_subtree(200, c200.begin(), c200.end(), sp);
        std::FILE *fp(tmpfile());
        REQUIRE(fp != nullptr);
        ks::string ks;
        fme.run_collapse(1000, ks, fp);
        // Only the two shared pairs in subtree 100 survive the bitmap filter;
        // a two-genome subtree yields no candidates at all.
        REQUIRE(fme.added_nodes().size() == 2);
        for(const auto &pair: fme.added_nodes()) REQUIRE(pair.first == 100);
        const auto emitted(parse_emitted(fp));
        std::fclose(fp);
        REQUIRE(emitted.size() == 2);
    }

    kh_destroy(p, taxmap);
}
