#ifndef BONSAI_PROFILE_STATS_H__
#define BONSAI_PROFILE_STATS_H__
// Statistics over sample profiles (bonsai profile diff) and the sample index
// for nearest-profile search (bonsai profile index / search).
//
// Differential abundance between two groups of samples, per clade:
//  - reads: Welch t on log proportions (log((reads + 0.5) / (classified + 0.5)): exactly 0
//    for a clade holding every read, whatever the depth);
//    p by permuting the group labels (all relabelings when there are at most
//    `nperm` of them, else nperm random ones; p = (1 + #|t*| >= |t|) / (1 +
//    nperm)), so no normality is assumed, only exchangeable samples under the
//    null. Benjamini-Hochberg q over the clades tested.
//  - k-mers: PERMANOVA (pseudo-F) on the samples' pairwise k-mer Jaccard
//    distances for that clade, permutation p as above: a change of the clade's
//    k-mer content (strain replacement, gene gain or loss) even when its
//    abundance does not move.
// And community-wide: PERMANOVA on tree EMD and on Bray-Curtis.
// Small groups limit the attainable p (two samples per group: 6 relabelings,
// p >= 1/3); the output prints the number of relabelings.
#include <random>
#include "profile.h"

namespace bns {

// Benjamini-Hochberg; NaN p-values (untested) stay NaN and are not counted.
inline std::vector<double> bh_qvalues(const std::vector<double> &p) {
    std::vector<size_t> idx;
    for(size_t i = 0; i < p.size(); ++i) if(!std::isnan(p[i])) idx.push_back(i);
    std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b) {return p[a] < p[b];});
    std::vector<double> q(p.size(), std::numeric_limits<double>::quiet_NaN());
    double run = 1.;
    for(size_t r = idx.size(); r-- > 0;) {
        run = std::min(run, p[idx[r]] * idx.size() / double(r + 1));
        q[idx[r]] = run;
    }
    return q;
}

inline double n_choose_k(unsigned n, unsigned k) {
    double r = 1.;
    for(unsigned i = 1; i <= k; ++i) r = r * (n - k + i) / i;
    return r;
}

// Permutation p for stat(labels) (labels[i] = 1: group b), larger = more
// extreme. Returns {p, relabelings used}.
template<typename Stat>
std::pair<double, size_t> permutation_p(const std::vector<uint8_t> &labels, Stat &&stat, size_t nperm, uint64_t seed) {
    const double obs = stat(labels);
    if(std::isnan(obs)) return {std::numeric_limits<double>::quiet_NaN(), 0};
    const double tol = 1e-12 * std::max(1., std::abs(obs));
    const unsigned nb = std::count(labels.begin(), labels.end(), 1);
    if(n_choose_k(labels.size(), nb) <= double(nperm)) {   // every relabeling, the observed one included
        std::vector<uint8_t> l(labels.size(), 0);
        std::fill(l.end() - nb, l.end(), 1);
        size_t ge = 0, n = 0;
        do {++n; ge += stat(l) >= obs - tol;} while(std::next_permutation(l.begin(), l.end()));
        return {double(ge) / n, n};
    }
    std::mt19937_64 rng(seed);
    std::vector<uint8_t> l = labels;
    size_t ge = 0;
    for(size_t i = 0; i < nperm; ++i) {
        std::shuffle(l.begin(), l.end(), rng);
        ge += stat(l) >= obs - tol;
    }
    return {(1. + ge) / (1. + nperm), nperm};
}

// |Welch t| of x between label groups (NaN if a group has < 2 values).
inline double welch_abs_t(const std::vector<double> &x, const std::vector<uint8_t> &l) {
    double s[2] = {0., 0.}, ss[2] = {0., 0.};
    size_t n[2] = {0, 0};
    for(size_t i = 0; i < x.size(); ++i) s[l[i]] += x[i], ss[l[i]] += x[i] * x[i], ++n[l[i]];
    if(n[0] < 2 || n[1] < 2) return std::numeric_limits<double>::quiet_NaN();
    double m[2], v[2];
    for(int g = 0; g < 2; ++g) m[g] = s[g] / n[g], v[g] = std::max(0., (ss[g] - n[g] * m[g] * m[g]) / (n[g] - 1));
    const double se = std::sqrt(v[0] / n[0] + v[1] / n[1]), d = std::abs(m[1] - m[0]);
    return se > 0. ? d / se: d > 0. ? std::numeric_limits<double>::infinity(): 0.;
}

// PERMANOVA pseudo-F for a distance matrix (row-major n x n) and two groups.
inline double permanova_f(const std::vector<double> &d, size_t n, const std::vector<uint8_t> &l) {
    size_t ng[2] = {0, 0};
    for(const auto x: l) ++ng[x];
    if(ng[0] < 1 || ng[1] < 1 || n < 3) return std::numeric_limits<double>::quiet_NaN();
    double sst = 0., ssw[2] = {0., 0.};
    for(size_t i = 0; i < n; ++i) for(size_t j = i + 1; j < n; ++j) {
        const double x = d[i * n + j] * d[i * n + j];
        sst += x;
        if(l[i] == l[j]) ssw[l[i]] += x;
    }
    sst /= n;
    const double sw = ssw[0] / ng[0] + ssw[1] / ng[1];
    if(sw <= 0.) return sst > sw ? std::numeric_limits<double>::infinity(): 0.;
    return (sst - sw) / (sw / double(n - 2));
}

// ---- sample index (profile index / search) ----
// Per sample: name, read totals, the taxon and clade read counts, and the
// sample's k-mer counter (union of its root clades' counters: every distinct
// classified k-mer). gzip of a header (MAGIC, VERSION, k, canonical, kind, p,
// fingerprint) then entries to end of file; `profile index -a` appends a
// gzip member of entries, which gzread reads as one stream.
struct IndexEntry {
    std::string name;
    uint64_t classified = 0, unclassified = 0;
    std::unordered_map<tax_t, uint64_t> taxon_reads, clade_reads;
    UniqCounter sample;
};

struct SampleIndex {
    static constexpr char MAGIC[8] = {'B', 'N', 'S', 'P', 'I', 'D', 'X', '\0'};
    static constexpr uint32_t VERSION = 1;
    Profile header;   // k, canonical, kind, p, fingerprint only
    std::vector<IndexEntry> entries;

    static void put_map(gzFile fp, const std::unordered_map<tax_t, uint64_t> &m) {
        std::vector<std::pair<tax_t, uint64_t>> v(m.begin(), m.end());
        std::sort(v.begin(), v.end());
        Profile::put(fp, uint64_t(v.size()));
        for(const auto &kv: v) Profile::put(fp, uint32_t(kv.first)), Profile::put(fp, kv.second);
    }
    static void get_map(gzFile fp, std::unordered_map<tax_t, uint64_t> &m) {
        for(uint64_t n = Profile::get<uint64_t>(fp); n--;) {
            const tax_t t = Profile::get<uint32_t>(fp);
            m[t] = Profile::get<uint64_t>(fp);
        }
    }
    static void write_header(gzFile fp, const Profile &h) {
        UniqCounter::gzwrite_all(fp, MAGIC, sizeof(MAGIC));
        Profile::put(fp, VERSION); Profile::put(fp, h.k);
        const uint8_t b[3] = {h.canonical, h.kind, h.p};
        UniqCounter::gzwrite_all(fp, b, sizeof(b));
        Profile::put(fp, h.fingerprint);
    }
    static void write_entry(gzFile fp, IndexEntry &e) {
        Profile::put(fp, uint32_t(e.name.size()));
        UniqCounter::gzwrite_all(fp, e.name.data(), e.name.size());
        Profile::put(fp, e.classified); Profile::put(fp, e.unclassified);
        put_map(fp, e.taxon_reads); put_map(fp, e.clade_reads);
        e.sample.write(fp);
    }
    // header only if entries_too is false (append checks)
    static SampleIndex read(const std::string &path, bool entries_too = true) {
        SampleIndex ret;
        gzFile fp = Profile::open(path, "rb");
        try {
            char magic[sizeof(MAGIC)];
            UniqCounter::gzread_all(fp, magic, sizeof(magic));
            if(std::memcmp(magic, MAGIC, sizeof(MAGIC))) RUNTIME_ERROR(path + " is not a bonsai sample index");
            if(Profile::get<uint32_t>(fp) != VERSION) RUNTIME_ERROR(path + ": unsupported index version");
            ret.header.k = Profile::get<uint32_t>(fp);
            uint8_t b[3];
            UniqCounter::gzread_all(fp, b, sizeof(b));
            ret.header.canonical = b[0]; ret.header.kind = b[1]; ret.header.p = b[2];
            ret.header.fingerprint = Profile::get<uint64_t>(fp);
            while(entries_too) {
                uint32_t len;
                const int r = gzread(fp, &len, sizeof(len));
                if(r == 0) break;
                if(r != int(sizeof(len))) RUNTIME_ERROR(path + ": truncated index");
                IndexEntry &e = ret.entries.emplace_back();
                e.name.resize(len);
                UniqCounter::gzread_all(fp, e.name.data(), len);
                e.classified = Profile::get<uint64_t>(fp); e.unclassified = Profile::get<uint64_t>(fp);
                get_map(fp, e.taxon_reads); get_map(fp, e.clade_reads);
                e.sample.read(fp);
            }
        } catch(...) {gzclose(fp); throw;}
        gzclose(fp);
        return ret;
    }
};

inline IndexEntry index_entry(ProfileSummary &&s) {
    IndexEntry e;
    e.name = std::move(s.name); e.classified = s.classified; e.unclassified = s.unclassified;
    e.taxon_reads = std::move(s.taxon_reads); e.clade_reads = std::move(s.clade_reads);
    e.sample = std::move(s.all);
    return e;
}

} // namespace bns

#endif /* BONSAI_PROFILE_STATS_H__ */
