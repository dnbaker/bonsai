/*
A/B bench: RangeMinHash<2048> vs BBMinHash (BBitMinHasher<uint64_t>, b = 16
retained bits per register) at equal register count (2048) and at equal
stored-sketch bytes (8192 registers x 16b = 16 KiB, the register content of
RMM-2048's 2048 x 64-bit minima).

Three parts, all over canonical k-mers hashed with the same sketch::WangHash
the sketchers apply internally (so exact Jaccard on the hashed sets is ground
truth; 64-bit collisions are negligible):

1. Real genomes (test/*.fa{,.gz}): sketch time per genome, final-jaccard
   time per pair, and |J_est - J_true| over all pairs.
2. H0 noise floor: iid-random 1500bp genomes (mirrors the
   sketch_candidates_nulls fixture). Reports the J-estimate distribution
   under the null and how many of the 496 pairs cross the tested
   thresholds (0.005 / 0.02).
3. Error-vs-J curve: 8kbp genomes with injected shared blocks at controlled
   J. Locates the matched-memory crossover where each estimator wins.

Timing is single-threaded; sketch construction reps take min/mean.
*/
#include "encoder.h"
#include "spacer.h"
#include "sketch/mh.h"
#include "sketch/bbmh.h"
#include "aesctr/wy.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <numeric>
#include <string>
#include <vector>

using namespace bns;

namespace {

using HSet = std::vector<u64>; // sorted unique WangHash64 of canonical k-mers

sketch::WangHash g_hf;

HSet exact_hashes(const std::string &path, const Spacer &sp) {
    std::vector<u64> raw;
    Encoder<score::Lex> enc(0, 0, sp, nullptr, true);
    enc.for_each([&](u64 v) {raw.push_back(g_hf(v));}, path.data());
    std::sort(raw.begin(), raw.end());
    raw.erase(std::unique(raw.begin(), raw.end()), raw.end());
    return raw;
}

HSet exact_hashes_seq(const std::string &seq, const Spacer &sp) {
    std::vector<u64> raw;
    Encoder<score::Lex> enc(0, 0, sp, nullptr, true);
    enc.assign(seq.data(), seq.size());
    enc.for_each_canon_unwindowed([&](u64 v) {raw.push_back(g_hf(v));});
    std::sort(raw.begin(), raw.end());
    raw.erase(std::unique(raw.begin(), raw.end()), raw.end());
    return raw;
}

double exact_jaccard(const HSet &a, const HSet &b) {
    size_t shared = 0;
    auto ia = a.begin(), ib = b.begin();
    while(ia != a.end() && ib != b.end()) {
        if(*ia < *ib) ++ia;
        else if(*ib < *ia) ++ib;
        else {++shared; ++ia; ++ib;}
    }
    return double(shared) / (a.size() + b.size() - shared);
}

inline std::string rand_dna(wy::WyRand<u64> &rng, size_t len) {
    static const char bases[] = "ACGT";
    std::string ret;
    ret.reserve(len);
    while(ret.size() < len) ret.push_back(bases[rng() & 3]);
    return ret;
}

template<typename F>
double time_ms(F &&f, unsigned reps, double &mean_out) {
    std::vector<double> ts;
    ts.reserve(reps);
    for(unsigned i(0); i < reps; ++i) {
        const auto t0(std::chrono::steady_clock::now());
        f();
        const auto t1(std::chrono::steady_clock::now());
        ts.push_back(std::chrono::duration<double, std::milli>(t1 - t0).count());
    }
    mean_out = std::accumulate(ts.begin(), ts.end(), 0.) / reps;
    return *std::min_element(ts.begin(), ts.end());
}

// ---- sketch adapters ----
struct RMM2048 {
    static constexpr const char *name = "RMM-2048x64b";
    static constexpr size_t live_nbytes = 2048 * 8;
    static constexpr size_t final_nbytes = 2048 * 8;
    sketch::RangeMinHash<u64> s;
    RMM2048(): s(2048) {}
    void addh(u64 v) {s.addh(v);}
    using final_type = sketch::RangeMinHash<u64>::final_type;
    final_type finalize() const {return s.finalize();}
};
struct BBMH8192 {
    static constexpr const char *name = "BBMH-8192regs b16";
    static constexpr size_t live_nbytes = size_t(8192) * 8;
    static constexpr size_t final_nbytes = size_t(16) * 8192 / 8;
    sketch::BBitMinHasher<u64> s;
    BBMH8192(): s(13, 16) {}
    void addh(u64 v) {s.addh(v);}
    using final_type = sketch::BBitMinHasher<u64>::final_type;
    final_type finalize() const {return s.finalize();}
};
struct BBMH2048 {
    static constexpr const char *name = "BBMH-2048regs b16";
    static constexpr size_t live_nbytes = size_t(2048) * 8;
    static constexpr size_t final_nbytes = size_t(16) * 2048 / 8;
    sketch::BBitMinHasher<u64> s;
    BBMH2048(): s(11, 16) {}
    void addh(u64 v) {s.addh(v);}
    using final_type = sketch::BBitMinHasher<u64>::final_type;
    final_type finalize() const {return s.finalize();}
};

template<typename Sketch>
void bench_real(const std::vector<std::string> &paths, const Spacer &sp,
                const std::vector<HSet> &truth) {
    std::fprintf(stderr, "== %s (live %zu B, final %zu B) ==\n",
                 Sketch::name, Sketch::live_nbytes, Sketch::final_nbytes);
    // one sketch per genome, untimed, for accuracy and jaccard timing
    std::vector<typename Sketch::final_type> fins;
    for(size_t i(0); i < paths.size(); ++i) {
        Sketch sk;
        Encoder<score::Lex> enc(0, 0, sp, nullptr, true);
        enc.for_each([&](u64 v) {sk.addh(v);}, paths[i].data());
        fins.push_back(sk.finalize());
    }
    double tot_min(0.), tot_mean(0.);
    volatile double sink0(0.);
    for(size_t i(0); i < paths.size(); ++i) {
        double mean;
        const double mn(time_ms([&]{
            Sketch sk;
            Encoder<score::Lex> enc(0, 0, sp, nullptr, true);
            enc.for_each([&](u64 v) {sk.addh(v);}, paths[i].data());
            sink0 += sk.finalize().jaccard_index(fins[i]);
        }, 5, mean));
        tot_min += mn; tot_mean += mean;
        std::fprintf(stderr, "  sketch %-40s %8.2f ms (mean %8.2f)\n",
                     paths[i].substr(paths[i].find_last_of('/') + 1).data(), mn, mean);
    }
    std::fprintf(stderr, "  sketch TOTAL min %8.2f ms  mean %8.2f ms\n", tot_min, tot_mean);
    const size_t np(paths.size() * (paths.size() - 1) / 2);
    volatile double sink(0.);
    double jmean;
    const unsigned jreps(20000 / np ? 20000 / np : 1);
    const double jmin(time_ms([&]{
        for(size_t i(0); i < paths.size(); ++i)
            for(size_t j(i + 1); j < paths.size(); ++j)
                sink += fins[i].jaccard_index(fins[j]);
    }, jreps, jmean));
    std::fprintf(stderr, "  jaccard %zu pairs x %u reps: %.3f ms/pass (%.0f ns/pair)\n",
                 np, jreps, jmin, jmin * 1e6 / np);
    double mae(0.), rmse(0.), mx(0.);
    std::fprintf(stderr, "  pair accuracy (J_true, J_est, |err|):\n");
    for(size_t i(0); i < paths.size(); ++i) {
        for(size_t j(i + 1); j < paths.size(); ++j) {
            const double jt(exact_jaccard(truth[i], truth[j]));
            const double je(fins[i].jaccard_index(fins[j]));
            const double ae(std::fabs(je - jt));
            mae += ae; rmse += ae * ae; mx = std::max(mx, ae);
            std::fprintf(stderr, "    %.6f  %.6f  %.6f\n", jt, je, ae);
        }
    }
    std::fprintf(stderr, "  MAE %.6f  RMSE %.6f  max %.6f\n\n",
                 mae / np, std::sqrt(rmse / np), mx);
}

template<typename Sketch>
void bench_null(const Spacer &sp) {
    wy::WyRand<u64> rng(99);
    std::vector<std::string> paths;
    std::vector<HSet> truth;
    for(unsigned g(0); g < 32; ++g) {
        const std::string seq(rand_dna(rng, 1500));
        truth.push_back(exact_hashes_seq(seq, sp));
        paths.push_back(std::move(seq));
    }
    std::vector<typename Sketch::final_type> fins;
    for(const auto &p: paths) {
        Sketch sk;
        Encoder<score::Lex> enc(0, 0, sp, nullptr, true);
        enc.assign(p.data(), p.size());
        enc.for_each_canon_unwindowed([&](u64 v) {sk.addh(v);});
        fins.push_back(sk.finalize());
    }
    const size_t n(32), np(n * (n - 1) / 2);
    double m(0.), sq(0.), mx(0.), mjt(0.);
    size_t ge5(0), ge20(0);
    for(size_t i(0); i < n; ++i)
        for(size_t j(i + 1); j < n; ++j) {
            const double je(fins[i].jaccard_index(fins[j]));
            const double jt(exact_jaccard(truth[i], truth[j]));
            m += je; sq += je * je; mx = std::max(mx, je); mjt += jt;
            ge5 += je >= 0.005; ge20 += je >= 0.02;
        }
    const double mean(m / np), sd(std::sqrt(sq / np - mean * mean));
    std::fprintf(stderr, "%-18s H0: mean %.5f sd %.5f max %.5f | mean J_true %.6f | pairs>=0.005: %zu/%zu  >=0.02: %zu/%zu\n",
                 Sketch::name, mean, sd, mx, mjt / np, ge5, np, ge20, np);
}

template<typename Sketch>
void bench_curve(const Spacer &sp, double jtarget, wy::WyRand<u64> &rng) {
    const size_t L(8000);
    std::vector<double> errs;
    for(unsigned rep(0); rep <12; ++rep) {
        const std::string a(rand_dna(rng, L));
        std::string b(rand_dna(rng, L));
        size_t todo(size_t(std::min(1., 2. * jtarget) * L) + 17);
        while(todo > 17) {
            const size_t blen(std::min(todo, size_t(40 + (rng() % 400))));
            const size_t apos(rng() % (L - blen)), bpos(rng() % (L - blen));
            b.replace(bpos, blen, a.substr(apos, blen));
            todo -= blen;
        }
        Sketch sa, sb;
        Encoder<score::Lex> enc(0, 0, sp, nullptr, true);
        enc.assign(a.data(), a.size());
        enc.for_each_canon_unwindowed([&](u64 v) {sa.addh(v);});
        enc.assign(b.data(), b.size());
        enc.for_each_canon_unwindowed([&](u64 v) {sb.addh(v);});
        const double jt(exact_jaccard(exact_hashes_seq(a, sp), exact_hashes_seq(b, sp)));
        if(jt <= 0) continue;
        errs.push_back(std::fabs(sa.finalize().jaccard_index(sb.finalize()) - jt) / jt);
    }
    const double mae(std::accumulate(errs.begin(), errs.end(), 0.) / errs.size());
    std::fprintf(stderr, "%-18s target J~%.4g: n %zu  mean |err|/J %.4f\n",
                 Sketch::name, jtarget, errs.size(), mae);
}

} // namespace

int main(int argc, char **argv) {
    Spacer sp(17, 17);
    const char *base(argc > 1 ? argv[1] : "test");
    std::vector<std::string> paths {
        std::string(base) + "/ec.k12.fa",
        std::string(base) + "/GCF_000302455.1_ASM30245v1_genomic.fna.gz",
        std::string(base) + "/GCF_000762265.1_ASM76226v1_genomic.fna.gz",
        std::string(base) + "/GCF_000953115.1_DSM1535_genomic.fna.gz",
        std::string(base) + "/phix.fa",
    };
    std::fprintf(stderr, "k=%u window=%u; canonical k-mers, WangHash64\n\n", sp.k_, sp.w_);
    std::vector<HSet> truth;
    for(const auto &p: paths) truth.push_back(exact_hashes(p, sp));
    std::fprintf(stderr, "Part 1: real genomes\n\n");
    bench_real<RMM2048>(paths, sp, truth);
    bench_real<BBMH8192>(paths, sp, truth);
    bench_real<BBMH2048>(paths, sp, truth);
    std::fprintf(stderr, "Part 2: H0 noise floor (32 iid-random 1500bp genomes, 496 pairs)\n\n");
    bench_null<RMM2048>(sp);
    bench_null<BBMH8192>(sp);
    bench_null<BBMH2048>(sp);
    std::fprintf(stderr, "\nPart 3: relative error vs J (8kbp genomes, injected sharing)\n\n");
    for(const double jt: {0.002, 0.008, 0.02, 0.05, 0.15, 0.5}) {
        wy::WyRand<u64> rng(u64(jt * 1e6) + 7);
        bench_curve<RMM2048>(sp, jt, rng);
        bench_curve<BBMH8192>(sp, jt, rng);
        bench_curve<BBMH2048>(sp, jt, rng);
        std::fprintf(stderr, "\n");
    }
    return 0;
}
