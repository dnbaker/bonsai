#include "test/catch.hpp"

#include "feature_min.h"
#include "util.h"

#include <map>
#include <random>
#include <set>
#include <string>
#include <vector>

using namespace bns;

namespace {

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

struct TmpDir {
    std::string dir;
    TmpDir(): dir([] {
        char tmpl[] = "/tmp/bns_freq_XXXXXX";
        const char *ret(mkdtemp(tmpl));
        REQUIRE(ret != nullptr);
        return std::string(ret);
    }()) {}
    ~TmpDir() {
        for(const auto &f: files) std::remove(f.data());
        std::remove(dir.data());
    }
    std::vector<std::string> files;
    std::string operator+(const char *suffix) const {return dir + "/" + suffix;}
    std::string add(const std::string &name, const std::string &content) {
        files.push_back(dir + "/" + name);
        return write_file(files.back(), content);
    }
};

// Like canon_kmers, but counts every occurrence (absolute frequency).
std::map<u64, u64> canon_kmer_counts(const std::string &seq, unsigned k) {
    auto val = [](char c) {return c == 'A' ? 0 : c == 'C' ? 1 : c == 'G' ? 2 : 3;};
    auto rc = [](char c) {return c == 'A' ? 'T' : c == 'T' ? 'A' : c == 'C' ? 'G' : 'C';};
    std::map<u64, u64> ret;
    for(size_t i(0); i + k <= seq.size(); ++i) {
        u64 f(0), rv(0);
        bool valid(true);
        for(unsigned j(0); j < k; ++j) {
            const char c(seq[i + j]);
            if(c != 'A' && c != 'C' && c != 'G' && c != 'T') {valid = false; break;}
            f = (f << 2) | val(c);
            rv = (rv << 2) | val(rc(seq[i + k - 1 - j]));
        }
        if(valid) ++ret[std::min(f, rv)];
    }
    return ret;
}

std::set<u64> canon_kmers(const std::string &seq, unsigned k) {
    std::set<u64> ret;
    for(const auto &pair: canon_kmer_counts(seq, k)) ret.insert(pair.first);
    return ret;
}

std::map<tax_t, tax_t> parent_pairs {{100, 1}, {101, 100}, {102, 100}, {103, 100}};

khash_t(p) *build_taxmap(const std::map<tax_t, tax_t> &pairs) {
    khash_t(p) *ret(kh_init(p));
    int khr;
    for(const auto &pair: pairs) {
        const auto ki(kh_put(p, ret, pair.first, &khr));
        kh_val(ret, ki) = pair.second;
    }
    {
        const auto ki(kh_put(p, ret, 1, &khr));
        kh_val(ret, ki) = 0;
    }
    return ret;
}

} // namespace

TEST_CASE("update_frequency_counter") {
    khash_t(p) *taxmap(build_taxmap(parent_pairs));
    // Genome 1: kmer A x3, kmer B x1. Genome 2: kmer A x2.
    khash_t(64) *g1(kh_init(64)), *g2(kh_init(64));
    auto bump = [](khash_t(64) *h, u64 key, u64 n) {
        int khr_;
        const auto it(kh_get(64, h, key));
        if(it == kh_end(h)) {
            const auto ki(kh_put(64, h, key, &khr_));
            kh_val(h, ki) = n;
        } else kh_val(h, it) += n;
    };
    bump(g1, 1000, 3); bump(g1, 2000, 1);
    bump(g2, 1000, 2);
    khash_t(64) *global(kh_init(64));
    update_frequency_counter(global, g1, taxmap, 101);
    update_frequency_counter(global, g2, taxmap, 102);
    REQUIRE(kh_size(global) == 2);
    const auto ia(kh_get(64, global, 1000)), ib(kh_get(64, global, 2000));
    REQUIRE(ia != kh_end(global));
    REQUIRE(ib != kh_end(global));
    REQUIRE(FMcount(kh_val(global, ia)) == 5); // absolute frequency: 3 + 2
    REQUIRE(FMtax(kh_val(global, ia)) == 100); // lca(101, 102)
    REQUIRE(FMcount(kh_val(global, ib)) == 1);
    REQUIRE(FMtax(kh_val(global, ib)) == 101);
    kh_destroy(64, global);
    kh_destroy(64, g1);
    kh_destroy(64, g2);
    kh_destroy(p, taxmap);
}

TEST_CASE("frequency_map") {
    TmpDir tmp;
    const std::string taxpath(tmp.add("tax.txt", "100|\t1\n101|\t100\n102|\t100\n103|\t100\n"));
    khash_t(p) *taxmap(build_parent_map(taxpath.data()));
    REQUIRE(taxmap != nullptr);

    const unsigned k(11);
    const std::string special("ACGTTGCAAGT");
    std::mt19937_64 rng(9001);
    auto genome = [&](unsigned n_special) {
        std::string ret(rand_dna(rng, 60));
        for(unsigned i(0); i < n_special; ++i) ret += special + rand_dna(rng, 17);
        return ret;
    };
    // g1: special x3 (tax 101), g2: special x2 (tax 102), g3: none (tax 103).
    const std::string g1seq(genome(3)), g2seq(genome(2)), g3seq(genome(0));
    tmp.add("g1.fa", ">g1\n" + g1seq + "\n");
    tmp.add("g2.fa", ">g2\n" + g2seq + "\n");
    tmp.add("g3.fa", ">g3\n" + g3seq + "\n");
    const std::string namepath(tmp.add("names.txt", "g1\t101\ng2\t102\ng3\t103\n"));

    // Brute-force ground truth: per-genome canonical k-mer counts.
    std::map<u64, u64> expected_count;
    std::map<u64, tax_t> expected_tax;
    const std::map<std::string, tax_t> seqs{{g1seq, 101}, {g2seq, 102}, {g3seq, 103}};
    for(const auto &pair: seqs) {
        for(const auto &kmer_pair: canon_kmer_counts(pair.first, k)) {
            const u64 key(kmer_pair.first);
            expected_count[key] += kmer_pair.second;
            const auto it(expected_tax.find(key));
            expected_tax[key] = it == expected_tax.end() ? pair.second
                                                           : lca(taxmap, it->second, pair.second);
        }
    }
    // The special k-mer must actually be frequent in g1/g2 and absent in g3.
    const auto special_kmers(canon_kmers(special, k));
    REQUIRE(special_kmers.size() == 1);
    const u64 special_key(*special_kmers.begin());
    REQUIRE(expected_count.at(special_key) == 5);
    REQUIRE(expected_tax.at(special_key) == 100);
    REQUIRE(expected_count.size() > 10);

    Spacer sp(k, k);
    const std::vector<std::string> paths{tmp.dir + "/g1.fa", tmp.dir + "/g2.fa", tmp.dir + "/g3.fa"};

    // Single-threaded and multi-threaded (counter reuse) builds must agree
    // with the brute-force counts entry for entry.
    for(int nt: {1, 2}) {
        khash_t(64) *map(frequency_map<score::Lex>(paths, taxmap, namepath.data(), sp, nt, true, 1 << 10));
        REQUIRE(map != nullptr);
        REQUIRE(kh_size(map) == expected_count.size());
        for(khiter_t ki(0); ki < kh_end(map); ++ki) {
            if(!kh_exist(map, ki)) continue;
            const auto it(expected_count.find(kh_key(map, ki)));
            INFO("k-mer present in frequency_map but not in ground truth");
            REQUIRE(it != expected_count.end());
            CHECK(FMcount(kh_val(map, ki)) == it->second);
            CHECK(FMtax(kh_val(map, ki)) == expected_tax.at(it->first));
        }
        kh_destroy(64, map);
    }
    kh_destroy(p, taxmap);
}
