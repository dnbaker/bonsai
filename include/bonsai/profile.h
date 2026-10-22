#ifndef BONSAI_PROFILE_H__
#define BONSAI_PROFILE_H__
// Sample profiles (classify -P <out.bprof>; bonsai profile report|merge|
// compare|table). A profile is one sample's classification result as
// sufficient statistics: per taxon its own reads and k-mer hits, its parent
// and rank, the database k-mers in its clade, and its distinct-k-mer counter
// (UniqCounter: sorted hashes, or raw unfinalized b-bit MinHash registers),
// written before the counters are folded up the tree. Every
// ancestor of a hit taxon is present, so a profile needs no taxonomy to be
// read, and clade totals, the classify -U report and per-clade sketches are
// all recomputed from it. Profiles are a commutative monoid: merging those of
// a sample's read subsets (lanes, chunks) gives byte for byte the profile of
// the whole, because counts add and counters are set unions (exact union,
// register max) in a canonical state (UniqCounter::canonicalize).
//
// File: gzip of MAGIC, u32 VERSION (2), u32 name length, name, u32 k, u8
// canonical, u8 counter kind (1: b-bit MinHash), u8 counter precision, u8
// has_ranks, u32 exact-set threshold (-e), u64 database
// fingerprint, u64 unclassified reads, u64 taxon count, then per taxon in
// ascending taxid: u32 taxid, u32 parent (0: none), u8 rank, u64 reads, u64
// k-mer hits, u64 clade database k-mers, counter. Little-endian host order.
// Version 1 (HyperLogLog counters) is rejected: regenerate those profiles.
#include "hashmap.h"
#include <cmath>
#include <cstring>
#include <map>
#include <set>
#include "uniq.h"

namespace bns {

inline const char *profile_rank_name(uint8_t r) {
    static const char *names[] = {"strain-level", "species", "genus", "family", "order", "class", "phylum", "kingdom", "domain"};
    return r <= 8 ? names[r]: "no rank";
}

struct ProfileTaxon {
    tax_t parent = 0;
    uint8_t rank = 255;
    uint64_t reads = 0, kmers = 0, db_kmers = 0;   // the taxon's own reads and hits; its clade's database k-mers
    UniqCounter uniq;
    ProfileTaxon() = default;
    ProfileTaxon(unsigned p, size_t exact_max): uniq(p, exact_max) {}
};

// A fingerprint of the database (k, canonicalization and its k-mer count per
// taxon): profiles are only comparable over the same database.
inline uint64_t db_fingerprint(const bns::umap<tax_t, uint64_t> &db_kmers, uint32_t k, bool canonical) {
    std::vector<std::pair<tax_t, uint64_t>> v(db_kmers.begin(), db_kmers.end());
    std::sort(v.begin(), v.end());
    uint64_t h = uniq_hash(uint64_t(k) << 1 | canonical);
    for(const auto &kv: v) h = uniq_hash(uniq_hash(h ^ kv.first) ^ kv.second);
    return h;
}

struct Profile {
    static constexpr char MAGIC[8] = {'B', 'N', 'S', 'P', 'R', 'O', 'F', '\0'};
    static constexpr uint32_t VERSION = 2;   // 2 adds the exact-set threshold; 1 (HLL counters) is rejected
    std::string name;
    uint32_t k = 0;
    uint8_t canonical = 1, kind = uint8_t(UniqSketch::BBMH), p = 0, has_ranks = 0;
    uint32_t exact_max = UniqCounter::EXACT_MAX;
    uint64_t fingerprint = 0, unclassified = 0;
    std::map<tax_t, ProfileTaxon> taxa;

    uint64_t classified() const {
        uint64_t n = 0;
        for(const auto &kv: taxa) n += kv.second.reads;
        return n;
    }
    template<typename T> static void put(gzFile fp, const T &v) {UniqCounter::gzwrite_all(fp, &v, sizeof(v));}
    template<typename T> static T get(gzFile fp) {T v; UniqCounter::gzread_all(fp, &v, sizeof(v)); return v;}
    static gzFile open(const std::string &path, const char *mode) {
        gzFile fp = gzopen(path.data(), mode);
        if(!fp) RUNTIME_ERROR("could not open profile " + path);
        return fp;
    }
    // Header; then the caller writes ntaxa records (write_taxon).
    void write_header(gzFile fp, uint64_t ntaxa) const {
        UniqCounter::gzwrite_all(fp, MAGIC, sizeof(MAGIC));
        put(fp, VERSION);
        put(fp, uint32_t(name.size()));
        UniqCounter::gzwrite_all(fp, name.data(), name.size());
        put(fp, k);
        const uint8_t b[4] = {canonical, kind, p, has_ranks};
        UniqCounter::gzwrite_all(fp, b, sizeof(b));
        put(fp, exact_max);
        put(fp, fingerprint); put(fp, unclassified); put(fp, ntaxa);
    }
    static void write_taxon(gzFile fp, tax_t t, tax_t parent, uint8_t rank, uint64_t reads, uint64_t kmers, uint64_t db_kmers, UniqCounter &u) {
        put(fp, uint32_t(t)); put(fp, uint32_t(parent)); put(fp, rank);
        put(fp, reads); put(fp, kmers); put(fp, db_kmers);
        u.write(fp);
    }
    void write(const std::string &path) {
        gzFile fp = open(path, "wb");
        write_header(fp, taxa.size());
        for(auto &kv: taxa)
            write_taxon(fp, kv.first, kv.second.parent, kv.second.rank, kv.second.reads, kv.second.kmers, kv.second.db_kmers, kv.second.uniq);
        if(gzclose(fp) != Z_OK) RUNTIME_ERROR("could not finish profile " + path);
    }
    static Profile read(const std::string &path) {
        Profile ret;
        gzFile fp = open(path, "rb");
        try {
            char magic[sizeof(MAGIC)];
            UniqCounter::gzread_all(fp, magic, sizeof(magic));
            if(std::memcmp(magic, MAGIC, sizeof(MAGIC))) RUNTIME_ERROR(path + " is not a bonsai profile");
            if(const uint32_t v = get<uint32_t>(fp); v != VERSION)
                RUNTIME_ERROR(path + ": unsupported profile version " + std::to_string(v) +
                              (v == 1 ? " (HyperLogLog-era profiles are no longer supported; regenerate with classify -P)": ""));
            ret.name.resize(get<uint32_t>(fp));
            UniqCounter::gzread_all(fp, ret.name.data(), ret.name.size());
            ret.k = get<uint32_t>(fp);
            uint8_t b[4];
            UniqCounter::gzread_all(fp, b, sizeof(b));
            ret.canonical = b[0]; ret.kind = b[1]; ret.p = b[2]; ret.has_ranks = b[3];
            if(UniqSketch(ret.kind) != UniqSketch::BBMH)
                RUNTIME_ERROR(path + ": a HyperLogLog profile (kind " + std::to_string(ret.kind) +
                              ") is no longer supported; regenerate with classify -P");
            ret.exact_max = get<uint32_t>(fp);
            ret.fingerprint = get<uint64_t>(fp); ret.unclassified = get<uint64_t>(fp);
            for(uint64_t n = get<uint64_t>(fp); n--;) {
                const tax_t t = get<uint32_t>(fp);
                ProfileTaxon &x = ret.taxa[t];
                x.parent = get<uint32_t>(fp); x.rank = get<uint8_t>(fp);
                x.reads = get<uint64_t>(fp); x.kmers = get<uint64_t>(fp); x.db_kmers = get<uint64_t>(fp);
                x.uniq = UniqCounter(ret.p, ret.exact_max);
                x.uniq.read(fp, ret.exact_max);
            }
        } catch(...) {gzclose(fp); throw;}
        gzclose(fp);
        return ret;
    }
    // Empty if comparable, else why not.
    std::string incompatibility(const Profile &o) const {
        if(k != o.k || canonical != o.canonical) return "different k or canonicalization";
        if(fingerprint != o.fingerprint) return "different databases";
        if(kind != o.kind || p != o.p) return "different distinct-k-mer counters (kind or -u precision)";
        if(exact_max != o.exact_max) return "different exact-set thresholds (-e)";
        return {};
    }
    // Union with another profile of the same database (reads of the same
    // sample, or a pooled group). Taxa present in both must agree on parent.
    void merge(Profile &&o) {
        if(const auto why = incompatibility(o); why.size()) RUNTIME_ERROR("cannot merge profiles " + name + " and " + o.name + ": " + why);
        unclassified += o.unclassified;
        has_ranks |= o.has_ranks;
        for(auto &kv: o.taxa) {
            auto it = taxa.find(kv.first);
            if(it == taxa.end()) {taxa.emplace(kv.first, std::move(kv.second)); continue;}
            ProfileTaxon &x = it->second, &y = kv.second;
            if(x.parent != y.parent) RUNTIME_ERROR("cannot merge profiles: taxon " + std::to_string(kv.first) + " has different parents");
            if(x.rank == 255) x.rank = y.rank;
            x.reads += y.reads; x.kmers += y.kmers;
            x.uniq.merge(y.uniq);
        }
    }
    // The parent table (sketch cht) over this profile's taxa.
    CompactTaxMap parents() const {
        CompactTaxMap m(taxa.size() + 1);
        for(const auto &kv: taxa) m.set_parent(kv.first, kv.second.parent);
        return m;
    }
    // The counters as a UniqTable (moved out of the profile), to fold clades
    // and write the classify -U report.
    std::unique_ptr<UniqTable> to_table() {
        auto ret = std::make_unique<UniqTable>(1u, unsigned(p), size_t(exact_max));
        if(unclassified) ret->at(0).reads = unclassified;
        for(auto &kv: taxa) {
            SharedTally &st = ret->at(kv.first);
            st.reads = kv.second.reads; st.kmers = kv.second.kmers;
            st.uniq = std::move(kv.second.uniq);
            st.sketched = st.uniq.sketched();
        }
        return ret;
    }
};

// classify -P: the live table's taxon-level counts and counters, streamed
// (no copy) before the table is folded. Taxa are those with tallies and all
// their ancestors; a parent missing from the taxonomy is written as 0.
template<typename RankOf>
void write_profile(const std::string &path, Profile &header, UniqTable &table, const CompactTaxMap &taxmap,
                   const bns::umap<tax_t, uint64_t> &db_kmers, RankOf &&rank_of) {
    if(table.merged()) RUNTIME_ERROR("write_profile: the table has already been folded (write the profile before -U)");
    table.flush();
    auto parent = [&](tax_t t) -> tax_t {
        const tax_t par = taxmap.parent(t);
        return par == tax_t(-1) || par == t ? 0: par;
    };
    std::set<tax_t> nodes;
    table.for_each([&](tax_t t, SharedTally &) {if(t) nodes.insert(t);});
    for(std::vector<tax_t> todo(nodes.begin(), nodes.end()); todo.size();) {
        const tax_t t = todo.back(); todo.pop_back();
        if(const tax_t par = parent(t); par && nodes.insert(par).second) todo.push_back(par);
    }
    const auto db = clade_db_kmers(db_kmers, taxmap, [&](tax_t a) {return nodes.count(a) != 0;});
    header.kind = uint8_t(UniqSketch::BBMH); header.p = uint8_t(table.p); header.exact_max = uint32_t(table.exact_max);
    if(const SharedTally *u = table.find(0)) header.unclassified = u->reads.load();
    gzFile fp = Profile::open(path, "wb");
    header.write_header(fp, nodes.size());
    UniqCounter empty(table.p, table.exact_max);
    for(const tax_t t: nodes) {
        SharedTally *st = table.find(t);
        const auto it = db.find(t);
        Profile::write_taxon(fp, t, parent(t), rank_of(t), st ? st->reads.load(): 0, st ? st->kmers.load(): 0,
                             it == db.end() ? 0: it->second, st ? st->uniq: empty);
    }
    if(gzclose(fp) != Z_OK) RUNTIME_ERROR("could not finish profile " + path);
}

// One sample's clade view: taxon and clade read counts, and the clade
// counters (folded up the tree) of the taxa `keep_counter` selects.
struct ProfileSummary {
    std::string name;
    uint64_t classified = 0, unclassified = 0;
    bns::umap<tax_t, uint64_t> taxon_reads, clade_reads;
    bns::umap<tax_t, double> clade_uniq;
    bns::umap<tax_t, UniqCounter> counters;
    UniqCounter all;   // every distinct classified k-mer (union of the root clades)
};

// Consumes the profile. keep_counter(taxon, clade reads, rank) -> keep that
// clade's counter in the summary.
template<typename Keep>
ProfileSummary summarize(Profile &&pr, Keep &&keep_counter) {
    ProfileSummary s;
    s.name = pr.name; s.unclassified = pr.unclassified; s.classified = pr.classified();
    bns::umap<tax_t, uint8_t> ranks;
    for(const auto &kv: pr.taxa) {
        ranks.emplace(kv.first, kv.second.rank);
        if(kv.second.reads) s.taxon_reads.emplace(kv.first, kv.second.reads);
    }
    const CompactTaxMap parents = pr.parents();
    auto table = pr.to_table();
    table->merge_up(parents);
    s.all = UniqCounter(pr.p, pr.exact_max);
    table->for_each([&](tax_t t, SharedTally &st) {
        if(!t) return;
        if(const tax_t par = parents.parent(t); par == 0 || par == t || par == tax_t(-1)) s.all.merge(st.uniq);
        s.clade_reads.emplace(t, st.clade_reads);
        s.clade_uniq.emplace(t, st.uniq.estimate());
        const auto r = ranks.find(t);
        if(keep_counter(t, st.clade_reads, r == ranks.end() ? uint8_t(255): r->second)) s.counters.emplace(t, std::move(st.uniq));
    });
    return s;
}

// Distinct k-mers shared by two counters of one precision: exact for two
// exact counters; otherwise the one-permutation estimate J = matching
// registers / registers set in either (registers keep 64 - p hash bits, so
// false matches are negligible). An exact counter meets a sketch as a sketch
// of its hashes.
struct CounterOverlap {double a = 0., b = 0., shared = 0., jaccard = 0.;};
inline CounterOverlap counter_overlap(UniqCounter &x, UniqCounter &y) {
    CounterOverlap r;
    r.a = x.estimate(); r.b = y.estimate();
    if(!x.sketched() && !y.sketched()) {
        const auto &u = x.exact(), &v = y.exact();
        size_t n = 0;
        for(size_t i = 0, j = 0; i < u.size() && j < v.size();) {
            if(u[i] < v[j]) ++i; else if(v[j] < u[i]) ++j; else ++n, ++i, ++j;
        }
        r.shared = double(n);
        const double un = r.a + r.b - r.shared;
        r.jaccard = un > 0. ? r.shared / un: 0.;
        return r;
    }
    auto regs = [](const UniqCounter &c, std::vector<uint64_t> &tmp) -> const std::vector<uint64_t> & {
        if(c.bbmh_sketch()) {
            tmp.assign(c.bbmh_sketch()->core().begin(), c.bbmh_sketch()->core().end());
            return tmp;
        }
        uniq_bbmh_t s(c.p());
        for(const auto h: c.exact()) s.add(h);
        tmp.assign(s.core().begin(), s.core().end());
        return tmp;
    };
    std::vector<uint64_t> ta, tb;
    const auto &ra = regs(x, ta), &rb = regs(y, tb);
    static constexpr uint64_t EMPTY = std::numeric_limits<uint64_t>::max();
    size_t match = 0, either = 0;
    for(size_t i = 0; i < ra.size(); ++i) {
        either += ra[i] != EMPTY || rb[i] != EMPTY;
        match += ra[i] == rb[i] && ra[i] != EMPTY;
    }
    r.jaccard = either ? double(match) / either: 0.;
    r.shared = r.jaccard / (1. + r.jaccard) * (r.a + r.b);
    return r;
}

// Abundance distances between two samples. Bray-Curtis and presence Jaccard
// over `counts` (taxon -> reads, already restricted to the compared taxa;
// presence: >= min_reads); tree_emd: earth mover's distance on the taxonomy
// with unit edges between the classified-read distributions, i.e. the sum over
// clades of |fraction_a - fraction_b| (unnormalized weighted UniFrac with
// unit branch lengths; taxa outside the taxonomy hang off a virtual root).
// jaccard: the presence Jaccard distance (1 - shared / union).
struct AbundanceDistance {double bray_curtis = 0., jaccard = 0., tree_emd = 0.; size_t shared = 0;};
inline AbundanceDistance abundance_distance(const bns::umap<tax_t, uint64_t> &a, const bns::umap<tax_t, uint64_t> &b,
                                            const bns::umap<tax_t, uint64_t> &clade_a, uint64_t na,
                                            const bns::umap<tax_t, uint64_t> &clade_b, uint64_t nb, uint64_t min_reads) {
    AbundanceDistance d;
    double sa = 0., sb = 0.;
    for(const auto &kv: a) sa += kv.second;
    for(const auto &kv: b) sb += kv.second;
    double mins = 0.;
    size_t uni = 0;
    for(const auto &kv: a) {
        const auto it = b.find(kv.first);
        const uint64_t vb = it == b.end() ? 0: it->second;
        if(sa > 0. && sb > 0.) mins += std::min(kv.second / sa, vb / sb);
        const bool pa = kv.second >= min_reads, pb = vb >= min_reads;
        d.shared += pa && pb; uni += pa || pb;
    }
    for(const auto &kv: b) if(!a.count(kv.first) && kv.second >= min_reads) ++uni;
    d.bray_curtis = sa > 0. && sb > 0. ? 1. - mins: double(sa > 0. || sb > 0.);
    d.jaccard = uni ? 1. - double(d.shared) / uni: 0.;
    auto frac = [](const bns::umap<tax_t, uint64_t> &m, tax_t t, uint64_t n) {
        const auto it = m.find(t);
        return it == m.end() || !n ? 0.: double(it->second) / n;
    };
    for(const auto &kv: clade_a) d.tree_emd += std::abs(frac(clade_a, kv.first, na) - frac(clade_b, kv.first, nb));
    for(const auto &kv: clade_b) if(!clade_a.count(kv.first)) d.tree_emd += frac(clade_b, kv.first, nb);
    return d;
}

} // namespace bns

#endif /* BONSAI_PROFILE_H__ */
