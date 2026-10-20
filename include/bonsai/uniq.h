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
// k-mers and a sketch beyond, so the many rarely hit taxa cost little: a
// HyperLogLog (sketch::hll_t, 2^p registers; -U) or a b-bit MinHash
// (sketch::BBitMinHasher, 2^p registers; -B), whose finalized per-clade
// sketches can be written and compared against reference sketches by
// containment or Jaccard. One shared counter per taxon (UniqTable): sketch
// adds are lock-free (atomic register min/max), exact counters take a
// per-taxon mutex, and each thread batches its exact adds (UniqBatch). After
// classification, merge_up folds every node into its parent, deepest first,
// so intermediate nodes hold their whole subtree. Without -U nothing is
// collected.
#include <algorithm>
#include <atomic>
#include <cstdio>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
#include "sketch/hll.h"
#include "sketch/bbmh.h"
#include "util.h"
#include "taxmap.h"

namespace bns {

inline uint64_t uniq_hash(uint64_t x) {   // splitmix64 finalizer
    x ^= x >> 30; x *= 0xbf58476d1ce4e5b9ULL;
    x ^= x >> 27; x *= 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}
inline uint64_t uniq_hash(u128 x) {return uniq_hash(uint64_t(x) ^ uniq_hash(uint64_t(x >> 64)));}

enum class UniqSketch: uint8_t {HLL, BBMH};

// Parent lookups for merge_up and the reports: a CompactTaxMap (sketch cht
// table), or a khash_t(p) for older callers; tax_t(-1) when absent.
inline tax_t tax_parent(const CompactTaxMap &m, tax_t t) {return m.parent(t);}
inline tax_t tax_parent(const khash_t(p) *m, tax_t t) {return get_parent(m, t);}
using uniq_bbmh_t = sketch::BBitMinHasher<uint64_t>;

class UniqCounter {
    std::vector<uint64_t> exact_;          // distinct hashes (unsorted pending tail past sorted_)
    size_t sorted_ = 0;
    std::unique_ptr<sketch::hll_t> hll_;
    std::unique_ptr<uniq_bbmh_t> bbmh_;
    unsigned p_;
    UniqSketch kind_;
    void compact() {
        std::sort(exact_.begin(), exact_.end());
        exact_.erase(std::unique(exact_.begin(), exact_.end()), exact_.end());
        sorted_ = exact_.size();
    }
    void promote() {
        if(kind_ == UniqSketch::HLL) {
            hll_.reset(new sketch::hll_t(p_));
            for(const auto h: exact_) hll_->add(h);
        } else {
            bbmh_.reset(new uniq_bbmh_t(p_));
            for(const auto h: exact_) bbmh_->add(h);
        }
        std::vector<uint64_t>().swap(exact_);
        sorted_ = 0;
    }
public:
    // EXACT_MAX hashes take 8 * 512 = 4 KiB, the size of a p = 12 HyperLogLog.
    // (b-bit MinHash estimates and finalized cardinalities are correct for
    // sparse sketches from sketch d45a3c3 on: linear counting while registers
    // are empty.)
    static constexpr size_t EXACT_MAX = 512;
    explicit UniqCounter(unsigned p = 12, UniqSketch kind = UniqSketch::HLL): p_(p), kind_(kind) {}
    void add(uint64_t h) {
        if(hll_) {hll_->add(h); return;}
        if(bbmh_) {bbmh_->add(h); return;}
        exact_.push_back(h);
        if(exact_.size() >= 2 * EXACT_MAX) {
            compact();
            if(exact_.size() > EXACT_MAX) promote();
        }
    }
    // Union. An exact counter absorbing a sketch takes the sketch's kind and
    // precision; two sketches must share both.
    void merge(UniqCounter &o) {
        if(o.sketched()) {
            if(!sketched()) p_ = o.p_, kind_ = o.kind_, promote();
            else if(p_ != o.p_ || kind_ != o.kind_)
                throw std::invalid_argument("UniqCounter::merge: sketches of different kind or precision");
            if(hll_) *hll_ += *o.hll_;
            else *bbmh_ += *o.bbmh_;
            return;
        }
        for(const auto h: o.exact_) add(h);
    }
    // Exact iff at most EXACT_MAX distinct hashes, a sketch beyond, so the
    // state (and every estimate) depends only on the set of hashes added, not
    // on the order adds arrived in or on merges (profiles rely on this).
    void canonicalize() {
        if(sketched()) return;
        compact();
        if(exact_.size() > EXACT_MAX) promote();
    }
    double estimate() {
        canonicalize();
        if(hll_) return hll_->report();
        if(bbmh_) return bbmh_->cardinality_estimate();
        return double(exact_.size());
    }
    bool sketched() const {return hll_ || bbmh_;}
    UniqCounter clone() const {
        UniqCounter r(p_, kind_);
        r.exact_ = exact_; r.sorted_ = sorted_;
        if(hll_) r.hll_.reset(new sketch::hll_t(*hll_));
        if(bbmh_) r.bbmh_.reset(new uniq_bbmh_t(*bbmh_));
        return r;
    }
    unsigned p() const {return p_;}
    // After canonicalize(): the sorted distinct hashes (empty once sketched),
    // or the sketch (null unless that kind).
    const std::vector<uint64_t> &exact() const {return exact_;}
    const sketch::hll_t *hll() const {return hll_.get();}
    const uniq_bbmh_t *bbmh_sketch() const {return bbmh_.get();}
    // Serialized (little-endian host order): kind and precision bytes, a state
    // byte (0 exact, 1 sketch), then either a u32 count and the sorted distinct
    // hashes, or the raw registers (2^p bytes for HyperLogLog, 2^p u64 for
    // b-bit MinHash, unfinalized), so a loaded counter merges like a live one.
    void write(gzFile fp) {
        canonicalize();
        const uint8_t hdr[3] = {uint8_t(kind_), uint8_t(p_), uint8_t(sketched())};
        gzwrite_all(fp, hdr, sizeof(hdr));
        if(hll_) gzwrite_all(fp, hll_->core().data(), hll_->core().size());
        else if(bbmh_) gzwrite_all(fp, bbmh_->core().data(), bbmh_->core().size() * sizeof(uint64_t));
        else {
            const uint32_t n = exact_.size();
            gzwrite_all(fp, &n, sizeof(n));
            gzwrite_all(fp, exact_.data(), n * sizeof(uint64_t));
        }
    }
    void read(gzFile fp) {
        release();
        uint8_t hdr[3];
        gzread_all(fp, hdr, sizeof(hdr));
        if(hdr[0] > uint8_t(UniqSketch::BBMH) || hdr[1] < 4 || hdr[1] > 24 || hdr[2] > 1)
            throw std::runtime_error("UniqCounter::read: malformed counter");
        kind_ = UniqSketch(hdr[0]); p_ = hdr[1];
        if(hdr[2] && kind_ == UniqSketch::HLL) {
            hll_.reset(new sketch::hll_t(p_));
            gzread_all(fp, hll_->mutable_core().data(), hll_->core().size());
            hll_->not_ready();
        } else if(hdr[2]) {
            // Registers hold hash bits below the top p; re-adding the register
            // index on top restores each one exactly.
            std::vector<uint64_t> regs(size_t(1) << p_);
            gzread_all(fp, regs.data(), regs.size() * sizeof(uint64_t));
            bbmh_.reset(new uniq_bbmh_t(p_));
            for(size_t i = 0; i < regs.size(); ++i)
                if(regs[i] != std::numeric_limits<uint64_t>::max()) bbmh_->add((uint64_t(i) << (64 - p_)) | regs[i]);
        } else {
            uint32_t n;
            gzread_all(fp, &n, sizeof(n));
            if(n > EXACT_MAX) throw std::runtime_error("UniqCounter::read: malformed counter");
            exact_.resize(n);
            gzread_all(fp, exact_.data(), n * sizeof(uint64_t));
            sorted_ = n;
        }
    }
    static void gzwrite_all(gzFile fp, const void *p, size_t n) {
        if(n && gzwrite(fp, p, unsigned(n)) != int(n)) throw std::runtime_error("UniqCounter: write failed");
    }
    static void gzread_all(gzFile fp, void *p, size_t n) {
        if(n && gzread(fp, p, unsigned(n)) != int(n)) throw std::runtime_error("UniqCounter: truncated input");
    }
    // A sketched counter only (see SharedTally): safe concurrently with
    // other sketched adds, since both sketches update registers by CAS.
    void add_concurrent(uint64_t h) {
        if(hll_) hll_->add_threadsafe(h);
        else bbmh_->add_threadsafe(h);
    }
    // Drop the contents (the report frees each taxon's counter once merged).
    void release() {hll_.reset(); bbmh_.reset(); std::vector<uint64_t>().swap(exact_); sorted_ = 0;}
    UniqSketch kind() const {return kind_;}
    // The counter as a b-bit MinHash sketch (an exact counter is sketched on
    // the fly), finalized at b bits per register.
    sketch::FinalBBitMinHash bbmh(unsigned b) const {
        if(bbmh_) return bbmh_->finalize(b);
        if(hll_) throw std::invalid_argument("UniqCounter::bbmh: a HyperLogLog counter cannot be converted");
        uniq_bbmh_t tmp(p_);
        for(const auto h: exact_) tmp.add(h);
        return tmp.finalize(b);
    }
};

// One tally per taxon, shared by all classification threads, so sketch
// memory does not grow with the thread count. Sketched counters take
// concurrent adds without a lock (HyperLogLog and b-bit MinHash registers are
// updated by compare-and-swap); only the exact phase (a vector of hashes, and
// its promotion to a sketch) holds the tally's mutex, and threads batch those
// adds (UniqBatch) so a flush takes each tally's lock once.
struct SharedTally {
    std::atomic<uint64_t> reads{0};     // reads called this taxon
    std::atomic<uint64_t> kmers{0};     // k-mer hits mapping to this taxon
    std::atomic<bool> sketched{false};  // set (release) once uniq is a sketch
    std::mutex m;                       // guards uniq while it is exact
    UniqCounter uniq;                   // the taxon's k-mers; its clade's after UniqTable::merge_up
    // Set by UniqTable::merge_up: the taxon's own totals, and its clade's
    // (the taxon and all its descendants).
    uint64_t taxon_reads = 0, taxon_kmers = 0, clade_reads = 0, clade_kmers = 0;
    double taxon_uniq = 0.;
    SharedTally(unsigned p, UniqSketch kind): uniq(p, kind) {}
    template<typename It>
    void add(It beg, It end) {
        if(sketched.load(std::memory_order_acquire)) {
            for(; beg != end; ++beg) uniq.add_concurrent(*beg);
            return;
        }
        std::lock_guard<std::mutex> lock(m);
        for(; beg != end; ++beg) uniq.add(*beg);   // CAS adds once sketched
        if(uniq.sketched()) sketched.store(true, std::memory_order_release);
    }
};

class UniqTable;

// A thread's working group: pointers to the tallies it has touched (the
// table's shard lock is taken once per thread per taxon) and pending adds
// to exact tallies, applied in groups of FLUSH_AT.
class UniqBatch {
    struct Local {SharedTally *t; uint64_t reads = 0, kmers = 0;};
    UniqTable *table_ = nullptr;
    std::unordered_map<tax_t, Local> cache_;
    std::vector<std::pair<SharedTally *, uint64_t>> pending_;
    Local *last_ = nullptr;     // hits come in runs of one taxon
    tax_t last_t_ = 0;
    Local &local(tax_t t) {
        if(last_ && t == last_t_) return *last_;
        return *(last_ = &local_slow(last_t_ = t));
    }
    Local &local_slow(tax_t t);
public:
    static constexpr size_t FLUSH_AT = 4096;
    void bind(UniqTable *t) {table_ = t;}
    void add_kmer(tax_t t, uint64_t h) {
        Local &l = local(t);
        ++l.kmers;
        if(l.t->sketched.load(std::memory_order_acquire)) l.t->uniq.add_concurrent(h);
        else {
            pending_.emplace_back(l.t, h);
            if(pending_.size() >= FLUSH_AT) flush_pending();
        }
    }
    void add_read(tax_t t) {++local(t).reads;}
    void flush_pending() {
        std::sort(pending_.begin(), pending_.end(), [](const auto &a, const auto &b) {return a.first < b.first;});
        std::vector<uint64_t> group;
        for(size_t i = 0; i < pending_.size();) {
            size_t j = i;
            group.clear();
            for(; j < pending_.size() && pending_[j].first == pending_[i].first; ++j) group.push_back(pending_[j].second);
            pending_[i].first->add(group.begin(), group.end());
            i = j;
        }
        pending_.clear();
    }
    // Everything into the shared tallies; the cache is dropped.
    void flush() {
        flush_pending();
        for(auto &kv: cache_) {
            kv.second.t->reads.fetch_add(kv.second.reads, std::memory_order_relaxed);
            kv.second.t->kmers.fetch_add(kv.second.kmers, std::memory_order_relaxed);
        }
        std::unordered_map<tax_t, Local>().swap(cache_);
        last_ = nullptr;
    }
};

class UniqTable {
    static constexpr unsigned NSHARDS = 64;
    struct Shard {std::mutex m; std::unordered_map<tax_t, SharedTally> taxa;};
    std::unique_ptr<Shard[]> shards_{new Shard[NSHARDS]};
    std::vector<UniqBatch> batches_;
    bool merged_ = false;
public:
    const unsigned p;
    const UniqSketch kind;
    UniqTable(unsigned nthreads, unsigned p = 12, UniqSketch kind = UniqSketch::HLL): batches_(std::max(1u, nthreads)), p(p), kind(kind) {
        for(auto &b: batches_) b.bind(this);
    }
    UniqTable(const UniqTable &) = delete;
    SharedTally &at(tax_t t) {   // stable: unordered_map nodes do not move
        Shard &s = shards_[uniq_hash(uint64_t(t)) % NSHARDS];
        std::lock_guard<std::mutex> lock(s.m);
        return s.taxa.try_emplace(t, p, kind).first->second;
    }
    SharedTally *find(tax_t t) {
        Shard &s = shards_[uniq_hash(uint64_t(t)) % NSHARDS];
        std::lock_guard<std::mutex> lock(s.m);
        auto it = s.taxa.find(t);
        return it == s.taxa.end() ? nullptr: &it->second;
    }
    UniqBatch &batch(unsigned tid) {return batches_.at(tid);}
    bool merged() const {return merged_;}
    // After classification (no thread is adding).
    void flush() {for(auto &b: batches_) b.flush();}
    // After classification: flush, then fold every node into its parent,
    // deepest first, so each intermediate node's counter and clade totals hold
    // its whole subtree; each tree edge is merged once, in place (no second
    // copy of the counters). Ancestors without tallies of their own are added.
    // Idempotent.
    template<typename Parents>
    void merge_up(const Parents &taxmap) {
        flush();
        if(merged_) return;
        auto parent = [&](tax_t t) -> tax_t {
            const tax_t par = tax_parent(taxmap, t);
            return par == tax_t(-1) || par == t ? 0: par;
        };
        std::vector<tax_t> nodes;
        for_each([&](tax_t t, SharedTally &) {if(t) nodes.push_back(t);});
        std::unordered_map<tax_t, unsigned> depth;
        for(size_t i = 0; i < nodes.size(); ++i) {   // nodes grows with missing ancestors
            unsigned d = 0;
            for(tax_t a = parent(nodes[i]); a && d < 4096; a = parent(a), ++d)
                if(!find(a)) {at(a); nodes.push_back(a);}
            depth[nodes[i]] = d;
        }
        std::sort(nodes.begin(), nodes.end(), [&](tax_t a, tax_t b) {
            return depth[a] != depth[b] ? depth[a] > depth[b]: a < b;
        });
        for(const tax_t t: nodes) {
            SharedTally &st = *find(t);
            st.taxon_reads = st.reads.load(); st.taxon_kmers = st.kmers.load();
            st.clade_reads += st.taxon_reads; st.clade_kmers += st.taxon_kmers;
            st.taxon_uniq = st.uniq.estimate();   // before any child is folded in
        }
        for(const tax_t t: nodes) {   // children before parents
            SharedTally &st = *find(t);
            if(const tax_t par = parent(t)) {
                SharedTally &pt = *find(par);
                pt.clade_reads += st.clade_reads; pt.clade_kmers += st.clade_kmers;
                pt.uniq.merge(st.uniq);
            }
        }
        merged_ = true;
    }
    template<typename F> void for_each(F &&f) {
        for(unsigned i = 0; i < NSHARDS; ++i) for(auto &kv: shards_[i].taxa) f(kv.first, kv.second);
    }
};

inline UniqBatch::Local &UniqBatch::local_slow(tax_t t) {
    auto it = cache_.find(t);
    if(it == cache_.end()) it = cache_.emplace(t, Local{&table_->at(t)}).first;
    return it->second;
}

// Report: one row per taxon with reads or k-mers in its clade, by clade reads
// (descending). db_clade: database k-mers per clade (absent or 0: coverage
// NA); rank_name: rank label per taxon. sketch_dir (optional, b-bit MinHash
// counters): each clade's read k-mer sketch is written to
// <sketch_dir>/<taxid>.bbmh, finalized at sketch_b bits per register.
template<typename Parents, typename RankName>
void write_uniq_rows(std::FILE *fp, UniqTable &table, const Parents &taxmap,
                     const std::unordered_map<tax_t, uint64_t> &db, RankName &&rank_name,
                     const std::string *sketch_dir = nullptr, unsigned sketch_b = 16) {
    table.merge_up(taxmap);
    uint64_t unclassified = 0, total_reads = 0;
    std::vector<std::pair<tax_t, SharedTally *>> rows;
    table.for_each([&](tax_t t, SharedTally &st) {
        if(!t) {unclassified = st.reads.load(); total_reads += unclassified;}
        else {total_reads += st.taxon_reads; rows.emplace_back(t, &st);}
    });
    std::sort(rows.begin(), rows.end(), [](const auto &a, const auto &b) {
        return a.second->clade_reads != b.second->clade_reads ? a.second->clade_reads > b.second->clade_reads: a.first < b.first;
    });
    std::fprintf(fp, "#taxid\tparent\trank\tpct_reads\tclade_reads\ttaxon_reads\tclade_kmers\tclade_uniq_kmers\ttaxon_uniq_kmers\tdup\tcov\n");
    std::fprintf(fp, "0\t0\tunclassified\t%.4f\t%llu\t%llu\t0\t0\t0\t0\tNA\n",
                 total_reads ? 100. * unclassified / total_reads: 0., (unsigned long long)unclassified, (unsigned long long)unclassified);
    for(auto &r: rows) {
        SharedTally &c = *r.second;
        const double u = c.uniq.estimate();
        const tax_t par = tax_parent(taxmap, r.first);
        std::fprintf(fp, "%u\t%u\t%s\t%.4f\t%llu\t%llu\t%llu\t%.0f\t%.0f\t%.3f\t", r.first, par == tax_t(-1) ? 0u: par,
                     rank_name(r.first).data(), total_reads ? 100. * c.clade_reads / total_reads: 0.,
                     (unsigned long long)c.clade_reads, (unsigned long long)c.taxon_reads, (unsigned long long)c.clade_kmers,
                     u, c.taxon_uniq, u > 0. ? double(c.clade_kmers) / u: 0.);
        const auto it = db.find(r.first);
        if(it != db.end() && it->second) std::fprintf(fp, "%.6f\n", std::min(1., u / double(it->second)));
        else std::fprintf(fp, "NA\n");
        if(sketch_dir && table.kind == UniqSketch::BBMH)
            c.uniq.bbmh(sketch_b).write((*sketch_dir + "/" + std::to_string(r.first) + ".bbmh").data());
    }
}

// Database k-mers per clade, for each taxon `keep(t)` accepts: direct counts
// summed up the tree.
template<typename Parents, typename Keep>
std::unordered_map<tax_t, uint64_t> clade_db_kmers(const std::unordered_map<tax_t, uint64_t> &db_kmers, const Parents &taxmap, Keep &&keep) {
    std::unordered_map<tax_t, uint64_t> db;
    for(const auto &kv: db_kmers)
        for(tax_t a = kv.first, guard = 0; a && a != tax_t(-1) && guard < 4096; ++guard) {
            if(keep(a)) db[a] += kv.second;
            const tax_t par = tax_parent(taxmap, a);
            if(par == a) break;
            a = par;
        }
    return db;
}

// db_kmers (optional): database k-mers per taxon (direct), for coverage.
template<typename Parents, typename RankName>
void write_uniq_report(std::FILE *fp, UniqTable &table, const Parents &taxmap,
                       const std::unordered_map<tax_t, uint64_t> *db_kmers, RankName &&rank_name,
                       const std::string *sketch_dir = nullptr, unsigned sketch_b = 16) {
    table.merge_up(taxmap);
    std::unordered_map<tax_t, uint64_t> db;
    if(db_kmers) db = clade_db_kmers(*db_kmers, taxmap, [&](tax_t a) {return table.find(a) != nullptr;});
    write_uniq_rows(fp, table, taxmap, db, rank_name, sketch_dir, sketch_b);
}

} // namespace bns

#endif /* BONSAI_UNIQ_H__ */
