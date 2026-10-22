#ifndef BONSAI_UNIQ_H__
#define BONSAI_UNIQ_H__
// KrakenUniq-style distinct k-mer counts per taxon (classify -U). Every
// k-mer hit of every read is added to the counter of the taxon it maps to in
// the database (its LCA), and every read to the taxon it is called; the report
// adds clade totals (counter unions up the tree), duplication (hits per
// distinct k-mer) and coverage (distinct k-mers seen / the clade's k-mers in
// the database). Few distinct k-mers with many reads, or low coverage, flag
// false-positive taxa that read counts alone do not.
// Counters are exact (a sorted vector of hashes) up to a threshold of distinct
// k-mers (EXACT_MAX by default; a constructor parameter) and a b-bit MinHash
// sketch beyond (sketch::BBitMinHasher, 2^p registers), whose finalized
// per-clade sketches can be written and compared against reference sketches by
// containment or Jaccard. One shared counter per taxon (UniqTable): sketch
// adds are lock-free (atomic register min), exact counters take a per-taxon
// mutex, and each thread batches its exact adds (UniqBatch). After
// classification, merge_up folds every node into its parent, deepest first,
// so intermediate nodes hold their whole subtree. Without -U/-P nothing is
// collected.
#include "hashmap.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <atomic>
#include <cstdio>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
#include "sketch/bbmh.h"
#include "util.h"
#include "taxmap.h"

namespace bns {

inline uint64_t uniq_hash(uint64_t x) {
    // sketch::WangHash: the sketch submodule's core 64-bit hash (its default
    // Hasher, so counter register placement keeps the same distribution
    // quality), and a BIJECTION - WangHash::inverse recovers the pre-image.
    // (Was splitmix64's finalizer: same mixing strength, but not invertible.)
    return sketch::WangHash::hash(x);
}
inline uint64_t uniq_hash(u128 x) {return uniq_hash(uint64_t(x) ^ uniq_hash(uint64_t(x >> 64)));}

// The distinct-k-mer counter kind byte in profiles and counter serializations.
// Only b-bit MinHash remains; the byte is kept so the format stays stable and
// old HyperLogLog state is recognized (and rejected) rather than misparsed.
enum class UniqSketch: uint8_t {HLL = 0, BBMH = 1};

using uniq_bbmh_t = sketch::BBitMinHasher<uint64_t>;

class UniqCounter {
    std::vector<uint64_t> exact_;          // distinct hashes (unsorted pending tail past sorted_)
    size_t sorted_ = 0;
    std::unique_ptr<uniq_bbmh_t> bbmh_;
    unsigned p_;
    size_t exact_max_;                     // promote to a sketch past this many distinct hashes
    void compact() {
        std::sort(exact_.begin(), exact_.end());
        exact_.erase(std::unique(exact_.begin(), exact_.end()), exact_.end());
        sorted_ = exact_.size();
    }
    void promote() {
        bbmh_.reset(new uniq_bbmh_t(p_));
        for(const auto h: exact_) bbmh_->add(h);
        std::vector<uint64_t>().swap(exact_);
        sorted_ = 0;
    }
public:
    // EXACT_MAX hashes take 8 * 512 = 4 KiB per counter. (b-bit MinHash
    // estimates and finalized cardinalities are correct for sparse sketches
    // from sketch d45a3c3 on: linear counting while registers are empty.)
    static constexpr size_t EXACT_MAX = 512;
    explicit UniqCounter(unsigned p = 12, size_t exact_max = EXACT_MAX): p_(p), exact_max_(exact_max) {
        if(exact_max_ > (1ull << 31)) throw std::invalid_argument("UniqCounter: exact_max too large");
    }
    size_t exact_max() const {return exact_max_;}
    void add(uint64_t h) {
        if(bbmh_) {bbmh_->add(h); return;}
        exact_.push_back(h);
        if(exact_.size() >= 2 * exact_max_) {
            compact();
            if(exact_.size() > exact_max_) promote();
        }
    }
    // Union. An exact counter absorbing a sketch takes the sketch's
    // precision; two sketches must share theirs.
    void merge(UniqCounter &o) {
        if(o.sketched()) {
            if(!sketched()) p_ = o.p_, promote();
            else if(p_ != o.p_)
                throw std::invalid_argument("UniqCounter::merge: sketches of different precision");
            *bbmh_ += *o.bbmh_;
            return;
        }
        for(const auto h: o.exact_) add(h);
    }
    // Exact iff at most exact_max_ distinct hashes, a sketch beyond, so the
    // state (and every estimate) depends only on the set of hashes added, not
    // on the order adds arrived in or on merges (profiles rely on this).
    void canonicalize() {
        if(sketched()) return;
        compact();
        if(exact_.size() > exact_max_) promote();
    }
    double estimate() {
        canonicalize();
        if(bbmh_) return bbmh_->cardinality_estimate();
        return double(exact_.size());
    }
    bool sketched() const {return bool(bbmh_);}
    UniqCounter clone() const {
        UniqCounter r(p_, exact_max_);
        r.exact_ = exact_; r.sorted_ = sorted_;
        if(bbmh_) r.bbmh_.reset(new uniq_bbmh_t(*bbmh_));
        return r;
    }
    unsigned p() const {return p_;}
    // After canonicalize(): the sorted distinct hashes (empty once sketched),
    // or the sketch.
    const std::vector<uint64_t> &exact() const {return exact_;}
    const uniq_bbmh_t *bbmh_sketch() const {return bbmh_.get();}
    // Serialized (little-endian host order): kind and precision bytes, a state
    // byte (0 exact, 1 sketch), then either a u32 count and the sorted distinct
    // hashes, or the raw registers (2^p u64, unfinalized), so a loaded counter
    // merges like a live one. HyperLogLog state (kind byte 0) is no longer
    // supported: the HLL counter is gone, so old profiles must be regenerated.
    void write(gzFile fp) {
        canonicalize();
        const uint8_t hdr[3] = {uint8_t(UniqSketch::BBMH), uint8_t(p_), uint8_t(sketched())};
        gzwrite_all(fp, hdr, sizeof(hdr));
        if(bbmh_) gzwrite_all(fp, bbmh_->core().data(), bbmh_->core().size() * sizeof(uint64_t));
        else {
            const uint32_t n = exact_.size();
            gzwrite_all(fp, &n, sizeof(n));
            gzwrite_all(fp, exact_.data(), n * sizeof(uint64_t));
        }
    }
    void read(gzFile fp, size_t exact_max = EXACT_MAX) {
        release();
        exact_max_ = exact_max;
        uint8_t hdr[3];
        gzread_all(fp, hdr, sizeof(hdr));
        if(hdr[0] == uint8_t(UniqSketch::HLL))
            throw std::runtime_error("UniqCounter::read: a HyperLogLog counter (kind 0) is no longer supported; regenerate the profile with classify -P");
        if(hdr[0] != uint8_t(UniqSketch::BBMH) || hdr[1] < 4 || hdr[1] > 24 || hdr[2] > 1)
            throw std::runtime_error("UniqCounter::read: malformed counter");
        p_ = hdr[1];
        if(hdr[2]) {
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
            if(n > exact_max_) throw std::runtime_error("UniqCounter::read: malformed counter");
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
        bbmh_->add_threadsafe(h);
    }
    // Drop the contents (the report frees each taxon's counter once merged).
    void release() {bbmh_.reset(); std::vector<uint64_t>().swap(exact_); sorted_ = 0;}
    // The counter as a b-bit MinHash sketch (an exact counter is sketched on
    // the fly), finalized at b bits per register.
    sketch::FinalBBitMinHash bbmh(unsigned b) const {
        if(bbmh_) return bbmh_->finalize(b);
        uniq_bbmh_t tmp(p_);
        for(const auto h: exact_) tmp.add(h);
        return tmp.finalize(b);
    }
};

// One tally per taxon, shared by all classification threads, so sketch
// memory does not grow with the thread count. Sketched counters take
// concurrent adds without a lock (b-bit MinHash registers are updated by
// compare-and-swap); only the exact phase (a vector of hashes, and its
// promotion to a sketch) holds the tally's mutex, and threads batch those
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
    SharedTally(unsigned p, size_t exact_max): uniq(p, exact_max) {}
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
// table's shard lock is taken once per thread per taxon). Pending adds to
// exact tallies accumulate per tally (hits arrive in runs of one taxon, so
// no global sort is needed: each Local flushes its own small buffer).
class UniqBatch {
    struct Local {
        SharedTally *t;
        uint64_t reads = 0, kmers = 0;
        std::vector<uint64_t> pend;      // hashes awaiting add() to this exact tally
        void flush() {
            if(pend.empty()) return;
            t->add(pend.begin(), pend.end());
            pend.clear();
        }
    };
    UniqTable *table_ = nullptr;
    bns::smap<tax_t, Local> cache_;
    Local *last_ = nullptr;     // hits come in runs of one taxon
    tax_t last_t_ = 0;
    Local &local(tax_t t) {
        if(last_ && t == last_t_) return *last_;
        return *(last_ = &local_slow(last_t_ = t));
    }
    Local &local_slow(tax_t t);
public:
    static constexpr size_t FLUSH_AT = 64;   // pending hashes per tally before a flush
    void bind(UniqTable *t) {table_ = t;}
    void add_kmer(tax_t t, uint64_t h) {add_kmers(t, h, 1);}
    // n hits of one k-mer (a lookup run of equal keys). The counter holds a
    // set, so one add covers the whole run; the hit counter takes all n.
    void add_kmers(tax_t t, uint64_t h, uint64_t n) {
        Local &l = local(t);
        l.kmers += n;
        if(l.t->sketched.load(std::memory_order_acquire)) {
            l.t->uniq.add_concurrent(h);
            return;
        }
        l.pend.push_back(h);
        if(l.pend.size() >= FLUSH_AT) l.flush();
    }
    void add_read(tax_t t) {++local(t).reads;}
    void flush_pending() {
        for(auto &kv: cache_) kv.second.flush();
    }
    // Everything into the shared tallies; the cache is dropped.
    void flush() {
        flush_pending();
        for(auto &kv: cache_) {
            kv.second.t->reads.fetch_add(kv.second.reads, std::memory_order_relaxed);
            kv.second.t->kmers.fetch_add(kv.second.kmers, std::memory_order_relaxed);
        }
        bns::smap<tax_t, Local>().swap(cache_);
        last_ = nullptr;
    }
};

class UniqTable {
    // One table, one mutex: the mutex is taken only when a thread first
    // sights a taxon (a per-thread UniqBatch caches the tally pointer);
    // steady-state adds are lock-free (atomic counters, CAS sketch adds,
    // per-tally batched exact adds). Sharding added 64 maps and 64 mutexes
    // for no measured win at -p 4 (see the unshard measurement in the
    // round-2 report), so it is gone.
    std::mutex m_;
    bns::smap<tax_t, SharedTally> taxa_;   // segmented_map: element addresses are cached by batches and reports
    std::vector<UniqBatch> batches_;
    bool merged_ = false;
public:
    const unsigned p;
    const size_t exact_max;
    UniqTable(unsigned nthreads, unsigned p = 12, size_t exact_max = UniqCounter::EXACT_MAX): batches_(std::max(1u, nthreads)), p(p), exact_max(exact_max) {
        for(auto &b: batches_) b.bind(this);
    }
    UniqTable(const UniqTable &) = delete;
    SharedTally &at(tax_t t) {   // stable: segmented_map elements never move
        std::lock_guard<std::mutex> lock(m_);
        return taxa_.try_emplace(t, p, exact_max).first->second;
    }
    SharedTally *find(tax_t t) {
        std::lock_guard<std::mutex> lock(m_);
        auto it = taxa_.find(t);
        return it == taxa_.end() ? nullptr: &it->second;
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
        bns::umap<tax_t, unsigned> depth;
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
        for(auto &kv: taxa_) f(kv.first, kv.second);
    }
};

inline UniqBatch::Local &UniqBatch::local_slow(tax_t t) {
    auto it = cache_.find(t);
    if(it == cache_.end()) it = cache_.emplace(t, Local{&table_->at(t)}).first;
    return it->second;
}

// Per-k-mer coverage of a clade (Poisson): with P of its k-mers present at
// coverage lambda, hits = P lambda and distinct = P (1 - e^-lambda), so the
// multiplicity hits / distinct = lambda / (1 - e^-lambda) gives lambda, and P =
// hits / lambda (the sylph-style effective coverage, from the counter and the hit
// count together). NaN when no duplicates are seen (multiplicity <= 1:
// lambda unidentifiable, coverage << 1).
inline double clade_lambda(double hits, double distinct) {
    if(distinct <= 0. || hits <= distinct * (1. + 1e-9)) return std::numeric_limits<double>::quiet_NaN();
    const double m = hits / distinct;
    double lo = 0., hi = std::max(1., 2. * m);
    for(int i = 0; i < 200; ++i) {
        const double mid = .5 * (lo + hi);
        (mid / -std::expm1(-mid) < m ? lo: hi) = mid;
    }
    return .5 * (lo + hi);
}

// Report: one row per taxon with reads or k-mers in its clade, by clade reads
// (descending). db_clade: database k-mers per clade (absent or 0: coverage
// NA); rank_name: rank label per taxon. sketch_dir (optional, b-bit MinHash
// counters): each clade's read k-mer sketch is written to
// <sketch_dir>/<taxid>.bbmh, finalized at sketch_b bits per register.
template<typename Parents, typename RankName>
void write_uniq_rows(std::FILE *fp, UniqTable &table, const Parents &taxmap,
                     const bns::umap<tax_t, uint64_t> &db, RankName &&rank_name,
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
    std::fprintf(fp, "#taxid\tparent\trank\tpct_reads\tclade_reads\ttaxon_reads\tclade_kmers\tclade_uniq_kmers\ttaxon_uniq_kmers\tdup\tcov"
                     "\tlambda\tpresent_kmers\tcontainment\n");
    std::fprintf(fp, "0\t0\tunclassified\t%.4f\t%llu\t%llu\t0\t0\t0\t0\tNA\tNA\tNA\tNA\n",
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
        const bool have_db = it != db.end() && it->second;
        if(have_db) std::fprintf(fp, "%.6f", std::min(1., u / double(it->second)));
        else std::fprintf(fp, "NA");
        // lambda, present k-mers, containment (present / database k-mers of the clade)
        const double lam = clade_lambda(double(c.clade_kmers), u);
        if(std::isnan(lam)) std::fprintf(fp, "\tNA\tNA\tNA\n");
        else {
            const double present = double(c.clade_kmers) / lam;
            if(have_db) std::fprintf(fp, "\t%.4f\t%.0f\t%.6f\n", lam, present, std::min(1., present / double(it->second)));
            else std::fprintf(fp, "\t%.4f\t%.0f\tNA\n", lam, present);
        }
        if(sketch_dir)
            c.uniq.bbmh(sketch_b).write((*sketch_dir + "/" + std::to_string(r.first) + ".bbmh").data());
    }
}

// Database k-mers per clade, for each taxon `keep(t)` accepts: direct counts
// summed up the tree.
template<typename Parents, typename Keep>
bns::umap<tax_t, uint64_t> clade_db_kmers(const bns::umap<tax_t, uint64_t> &db_kmers, const Parents &taxmap, Keep &&keep) {
    bns::umap<tax_t, uint64_t> db;
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
                       const bns::umap<tax_t, uint64_t> *db_kmers, RankName &&rank_name,
                       const std::string *sketch_dir = nullptr, unsigned sketch_b = 16) {
    table.merge_up(taxmap);
    bns::umap<tax_t, uint64_t> db;
    if(db_kmers) db = clade_db_kmers(*db_kmers, taxmap, [&](tax_t a) {return table.find(a) != nullptr;});
    write_uniq_rows(fp, table, taxmap, db, rank_name, sketch_dir, sketch_b);
}

} // namespace bns

#endif /* BONSAI_UNIQ_H__ */
