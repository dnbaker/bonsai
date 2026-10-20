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

// Regularized incomplete beta I_x(a, b) (continued fraction, Lentz): binomial
// tails and Clopper-Pearson intervals for taxdiff.
inline double incbeta(double a, double b, double x) {
    if(x <= 0.) return 0.;
    if(x >= 1.) return 1.;
    if(x > (a + 1.) / (a + b + 2.)) return 1. - incbeta(b, a, 1. - x);
    const double lbeta = std::lgamma(a) + std::lgamma(b) - std::lgamma(a + b);
    const double front = std::exp(a * std::log(x) + b * std::log1p(-x) - lbeta) / a;
    double f = 1., c = 1., d = 0.;
    for(int i = 0; i <= 400; ++i) {
        const int m = i / 2;
        double num;
        if(i == 0) num = 1.;
        else if(i % 2 == 0) num = m * (b - m) * x / ((a + 2. * m - 1.) * (a + 2. * m));
        else num = -((a + m) * (a + b + m) * x) / ((a + 2. * m) * (a + 2. * m + 1.));
        d = 1. + num * d; if(std::abs(d) < 1e-300) d = 1e-300; d = 1. / d;
        c = 1. + num / c; if(std::abs(c) < 1e-300) c = 1e-300;
        const double cd = c * d;
        f *= cd;
        if(std::abs(1. - cd) < 1e-14) break;
    }
    return front * (f - 1.);
}
// x with incbeta(a, b, x) = q (bisection).
inline double incbeta_inv(double a, double b, double q) {
    double lo = 0., hi = 1.;
    for(int i = 0; i < 100; ++i) {
        const double mid = .5 * (lo + hi);
        (incbeta(a, b, mid) < q ? lo: hi) = mid;
    }
    return .5 * (lo + hi);
}
// Clopper-Pearson interval for a binomial proportion, k of n.
inline std::pair<double, double> clopper_pearson(double k, double n, double level = 0.95) {
    const double a = (1. - level) / 2.;
    return {k <= 0. ? 0.: incbeta_inv(k, n - k + 1., a), k >= n ? 1.: incbeta_inv(k + 1., n - k, 1. - a)};
}
// G statistic (2 x log-likelihood ratio) of k of n against expected proportion p0;
// chi-square 1 df p-value.
inline double binomial_g(double k, double n, double p0) {
    auto term = [](double o, double e) {return o > 0. ? o * std::log(o / e): 0.;};
    return n > 0. ? 2. * (term(k, n * p0) + term(n - k, n * (1. - p0))): 0.;
}
inline double chisq1_p(double g) {return std::erfc(std::sqrt(std::max(0., g) / 2.));}

// log beta-binomial pmf of k of n with mean p and intra-class correlation rho
// (a = p (1/rho - 1), b = (1 - p)(1/rho - 1)); the binomial when rho -> 0.
inline double log_betabinom(double k, double n, double p, double rho) {
    const double lc = std::lgamma(n + 1.) - std::lgamma(k + 1.) - std::lgamma(n - k + 1.);
    p = std::min(std::max(p, 1e-300), 1. - 1e-16);
    if(rho <= 1e-12) return lc + k * std::log(p) + (n - k) * std::log1p(-p);
    const double phi = 1. / rho - 1., a = p * phi, b = (1. - p) * phi;
    auto lbeta = [](double x, double y) {return std::lgamma(x) + std::lgamma(y) - std::lgamma(x + y);};
    return lc + lbeta(k + a, n - k + b) - lbeta(a, b);
}
// log e-value for k of n against null share p0 (beta-binomial, rho): the Bayes
// factor of a uniform mixture over the share against the null. With the null
// simple its expectation under the null is 1 (rho plugged in: approximate).
// The mixture integral is taken on a fine grid around the observed share.
inline double log_evalue_share(double k, double n, double p0, double rho) {
    if(n <= 0.) return 0.;
    const double lnull = log_betabinom(k, n, p0, rho);
    double lalt;
    if(rho <= 1e-12) lalt = -std::log(n + 1.);   // integral of the binomial over a uniform share
    else {
        const double c = (k + .5) / (n + 1.), sd = std::sqrt(c * (1. - c) * (rho + 1. / (n + 1.)));
        const double lo = std::max(1e-12, c - 12. * sd), hi = std::min(1. - 1e-12, c + 12. * sd);
        const int m = 800;
        const double h = (hi - lo) / m;
        double mx = -INFINITY;
        std::vector<double> v(m + 1);
        for(int i = 0; i <= m; ++i) mx = std::max(mx, v[i] = log_betabinom(k, n, lo + i * h, rho));
        double acc = 0.;
        for(int i = 0; i <= m; ++i) acc += (i == 0 || i == m ? .5: 1.) * std::exp(v[i] - mx);
        lalt = mx + std::log(acc * h);
    }
    return lalt - lnull;
}
// log e-value against an interval null: the share lies in [plo, phi] (e.g. a
// rate ratio within 2^-delta .. 2^delta). Null likelihood: its supremum over the
// band (beta-binomial, rho), which bounds every null in it, so the ratio stays
// a valid e-value; alternative: the uniform-share mixture of log_evalue_share.
inline double log_evalue_share_band(double k, double n, double plo, double phi, double rho) {
    if(n <= 0.) return 0.;
    const double lalt = log_evalue_share(k, n, 0.5, rho) + log_betabinom(k, n, 0.5, rho);   // the mixture's log marginal
    const double phat = std::min(std::max(k / n, plo), phi);
    return lalt - log_betabinom(k, n, phat, rho);
}

// e-BH (Wang & Ramdas 2022): the k largest e-values with e_(k) >= m / (alpha k)
// are rejected; valid FDR control under arbitrary dependence. Returns the
// threshold on e (infinity when nothing is rejected).
inline double ebh_threshold(std::vector<double> e, double alpha) {
    std::sort(e.begin(), e.end(), std::greater<double>());
    const double m = e.size();
    double t = INFINITY;
    for(size_t k = 1; k <= e.size(); ++k) if(e[k - 1] >= m / (alpha * k)) t = m / (alpha * k);
    return t;
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
