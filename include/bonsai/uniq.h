#ifndef BONSAI_UNIQ_H__
#define BONSAI_UNIQ_H__
// KrakenUniq-style distinct k-mer counts per taxon (classify -U). Every
// k-mer hit of every read is added to the counter of the taxon it maps to in
// the database (its LCA), and every read to the taxon it is called; the report
// adds clade totals (counter unions up the tree), duplication (hits per
// distinct k-mer) and coverage (distinct k-mers seen / the clade's k-mers in
// the database). Few distinct k-mers with many reads, or low coverage, flag
// false-positive taxa that read counts alone do not.
// Counters are exact (a sorted vector of hashes) up to EXACT_MAX distinct
// k-mers and a HyperLogLog (sketch::hll_t, 2^p registers) beyond, so the many
// rarely hit taxa cost little. One accumulator per classification thread,
// merged when the report is written; without -U nothing is collected.
#include <algorithm>
#include <cstdio>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>
#include "sketch/hll.h"
#include "util.h"

namespace bns {

inline uint64_t uniq_hash(uint64_t x) {   // splitmix64 finalizer
    x ^= x >> 30; x *= 0xbf58476d1ce4e5b9ULL;
    x ^= x >> 27; x *= 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}
inline uint64_t uniq_hash(u128 x) {return uniq_hash(uint64_t(x) ^ uniq_hash(uint64_t(x >> 64)));}

class UniqCounter {
    std::vector<uint64_t> exact_;          // distinct hashes (unsorted pending tail past sorted_)
    size_t sorted_ = 0;
    std::unique_ptr<sketch::hll_t> hll_;
    unsigned p_;
    void compact() {
        std::sort(exact_.begin(), exact_.end());
        exact_.erase(std::unique(exact_.begin(), exact_.end()), exact_.end());
        sorted_ = exact_.size();
    }
    void promote() {
        hll_.reset(new sketch::hll_t(p_));
        for(const auto h: exact_) hll_->add(h);
        std::vector<uint64_t>().swap(exact_);
        sorted_ = 0;
    }
public:
    // EXACT_MAX hashes take 8 * 512 = 4 KiB, the size of a p = 12 sketch
    static constexpr size_t EXACT_MAX = 512;
    explicit UniqCounter(unsigned p = 12): p_(p) {}
    void add(uint64_t h) {
        if(hll_) {hll_->add(h); return;}
        exact_.push_back(h);
        if(exact_.size() >= 2 * EXACT_MAX) {
            compact();
            if(exact_.size() > EXACT_MAX) promote();
        }
    }
    void merge(UniqCounter &o) {
        if(o.hll_) {
            if(!hll_) promote();
            *hll_ += *o.hll_;
            return;
        }
        for(const auto h: o.exact_) add(h);
    }
    double estimate() {
        if(hll_) return hll_->report();
        compact();
        return double(exact_.size());
    }
    bool sketched() const {return bool(hll_);}
};

struct TaxonTally {
    uint64_t reads = 0;      // reads called this taxon
    uint64_t kmers = 0;      // k-mer hits mapping to this taxon
    UniqCounter uniq;
    explicit TaxonTally(unsigned p = 12): uniq(p) {}
};

struct UniqAccumulator {
    std::unordered_map<tax_t, TaxonTally> taxa;
    unsigned p = 12;
    TaxonTally &at(tax_t t) {
        auto it = taxa.find(t);
        if(it == taxa.end()) it = taxa.emplace(t, TaxonTally(p)).first;
        return it->second;
    }
    void merge(UniqAccumulator &o) {
        for(auto &kv: o.taxa) {
            auto &t = at(kv.first);
            t.reads += kv.second.reads; t.kmers += kv.second.kmers;
            t.uniq.merge(kv.second.uniq);
        }
    }
};

// Report: one row per taxon with reads or k-mers in its clade, by clade reads
// (descending). db_kmers (optional): database k-mers per taxon (direct), for
// coverage; rank_name (optional): rank label per taxon.
template<typename RankName>
void write_uniq_report(std::FILE *fp, std::vector<UniqAccumulator> &accs, const khash_t(p) *taxmap,
                       const std::unordered_map<tax_t, uint64_t> *db_kmers, RankName &&rank_name) {
    UniqAccumulator all;
    if(!accs.empty()) all.p = accs.front().p;
    for(auto &a: accs) all.merge(a);
    struct Clade {uint64_t reads = 0, taxon_reads = 0, kmers = 0, db = 0; UniqCounter uniq; double taxon_uniq = 0.;
                  explicit Clade(unsigned p): uniq(p) {}};
    std::unordered_map<tax_t, Clade> clades;
    auto clade = [&](tax_t t) -> Clade & {
        auto it = clades.find(t);
        if(it == clades.end()) it = clades.emplace(t, Clade(all.p)).first;
        return it->second;
    };
    auto ancestors = [&](tax_t t, auto &&f) {
        for(unsigned guard = 0; t && t != tax_t(-1) && guard < 4096; ++guard) {
            f(t);
            t = get_parent(taxmap, t);
        }
    };
    uint64_t unclassified = 0, total_reads = 0;
    for(auto &kv: all.taxa) {
        total_reads += kv.second.reads;
        if(!kv.first) {unclassified += kv.second.reads; continue;}
        clade(kv.first).taxon_reads += kv.second.reads;
        clade(kv.first).taxon_uniq = kv.second.uniq.estimate();
        ancestors(kv.first, [&](tax_t a) {
            auto &c = clade(a);
            c.reads += kv.second.reads;
            c.kmers += kv.second.kmers;
            c.uniq.merge(kv.second.uniq);
        });
    }
    if(db_kmers)
        for(const auto &kv: *db_kmers)
            ancestors(kv.first, [&](tax_t a) {
                auto it = clades.find(a);
                if(it != clades.end()) it->second.db += kv.second;
            });
    std::vector<std::pair<tax_t, Clade *>> rows;
    for(auto &kv: clades) rows.emplace_back(kv.first, &kv.second);
    std::sort(rows.begin(), rows.end(), [](const auto &a, const auto &b) {
        return a.second->reads != b.second->reads ? a.second->reads > b.second->reads: a.first < b.first;
    });
    std::fprintf(fp, "#taxid\tparent\trank\tpct_reads\tclade_reads\ttaxon_reads\tclade_kmers\tclade_uniq_kmers\ttaxon_uniq_kmers\tdup\tcov\n");
    std::fprintf(fp, "0\t0\tunclassified\t%.4f\t%llu\t%llu\t0\t0\t0\t0\tNA\n",
                 total_reads ? 100. * unclassified / total_reads: 0., (unsigned long long)unclassified, (unsigned long long)unclassified);
    for(auto &r: rows) {
        Clade &c = *r.second;
        const double u = c.uniq.estimate();
        const tax_t par = get_parent(taxmap, r.first);
        std::fprintf(fp, "%u\t%u\t%s\t%.4f\t%llu\t%llu\t%llu\t%.0f\t%.0f\t%.3f\t", r.first, par == tax_t(-1) ? 0u: par,
                     rank_name(r.first).data(), total_reads ? 100. * c.reads / total_reads: 0.,
                     (unsigned long long)c.reads, (unsigned long long)c.taxon_reads, (unsigned long long)c.kmers,
                     u, c.taxon_uniq, u > 0. ? double(c.kmers) / u: 0.);
        if(db_kmers && c.db) std::fprintf(fp, "%.6f\n", std::min(1., u / double(c.db)));
        else std::fprintf(fp, "NA\n");
    }
}

} // namespace bns

#endif /* BONSAI_UNIQ_H__ */
