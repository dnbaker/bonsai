#include "test/catch.hpp"
#include "sketch/hll.h"
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <string>
#include <unordered_set>
#include <thread>
#include "encoder.h"
#include "omp.h"

using namespace bns;

// This file tests pure hll functionality
// and cardinal estimation, including
// a test for multi-threaded use.

int within_bounds(hll::hll_t &hll, double val) {
    return std::abs(hll.report() - val) <= hll.est_err();
}

std::vector<std::string> paths {
    "test/GCF_000302455.1_ASM30245v1_genomic.fna.gz",
    "test/GCF_000762265.1_ASM76226v1_genomic.fna.gz",
    "test/GCF_000953115.1_DSM1535_genomic.fna.gz",
    "test/GCF_001723155.1_ASM172315v1_genomic.fna.gz"
};

TEST_CASE("hll") {
    std::vector<std::pair<uint64_t, uint64_t>> pairs;
    {
        size_t sz_start /* = 9 */ = 10;
        size_t sz_end /* = 20 */ = 16;
        size_t nel_start = /* = 14 */ 18;
        size_t nel_stop = /* = 22 */ 21;
        for(size_t ss(sz_start); ss < sz_end; ++ss)
            for(size_t nel(nel_start); nel <= nel_stop; ++nel)
                pairs.emplace_back(std::pair<uint64_t, uint64_t>(ss, nel));
    }
    double diffsumsum = 0.;
    std::atomic<int> numpass; numpass.store(0);
    std::mt19937_64 mt;
    const size_t niter = 100;
    // const double div = niter;
    //std::fprintf(stdout, "#Structure Type\tSketch size\tNumber of elements\tDifference\tRMSE\tAbsDiffMean\tNumber underestimated\tNumber overestimated\tNumber of sketches with difference above expected\n");
    const auto npairs = pairs.size();
    omp_set_num_threads(std::thread::hardware_concurrency());
    std::set<std::pair<uint64_t, uint64_t>> gset;
    for(size_t ind = 0; ind < npairs; ++ind) {
        const auto pair(pairs[ind]);
        double diffsum = 0., absdiffsum = 0.;
        std::vector<double> diffs(niter);
        size_t numlessmore[] {0, 0};
        double sumlessmore[] {0, 0};
        int localnp = 0;
        OMP_PRAGMA("omp parallel for")
        for(size_t j = 0; j < niter; ++j) {
            std::set<std::pair<uint64_t, uint64_t>> lset;
            RNGType gen(mt() + j * (j + 1)); // I'm okay with a race condition, because I'm salting with the value of j
            const size_t val(1ull << pair.second), hsz(pair.first);
            hll::hll_t hll(hsz);
            for(size_t i(0); i < val; ++i) {
                hll.addh(gen());
            }
            hll.sum();
            if(std::abs(hll.report() - val) > hll.est_err()) {
                lset.insert(pair);
                //fprintf(stderr, "Warning: Above expected variance for %u, %u.\n", (unsigned)pair.first, (unsigned)pair.second);
            }
            auto diff = hll.report() - val; // For maximally mixed metaphors
            diffs[j] = diff;
            {
                OMP_PRAGMA("omp critical")
                diffsum += diff;
                ++numlessmore[diff > 0];
                sumlessmore[diff > 0] += std::abs(diff);
                absdiffsum += std::abs(diff);
                diff *= diff;
                localnp += (std::abs(hll.report() - val) <= hll.est_err());
                gset.insert(lset.begin(), lset.end());
            }
        }

        // std::fprintf(stderr, "Failed to meet expectations %zu of %zu times\n", gset.size(), niter);
        // auto sqvar = std::sqrt(std::accumulate(std::begin(diffs), std::end(diffs), 0., [div](auto x, auto y) {return x + y / div * y;}));
        {
//            std::fprintf(stdout, "HLL\t%u\t%u\t%lf\t%lf\t%lf\t%zu\t%zu\t%i\t%f\t%f\n",
//                         (unsigned)pair.first, (unsigned)pair.second, diffsum / div, sqvar, absdiffsum / div, numlessmore[0], numlessmore[1],
//                         (int)niter - localnp, sumlessmore[0] / div, sumlessmore[1] / div);
            diffsumsum += diffsum;
        }
        numpass += localnp;
    }
    gset.clear();
    for(size_t ind = 0; ind < npairs; ++ind) {
        const auto pair(pairs[ind]);
        double diffsum = 0., absdiffsum = 0.;
        std::vector<double> diffs(niter);
        size_t numlessmore[] {0, 0};
        double sumlessmore[] {0, 0};
        int localnp = 0;
        OMP_PRAGMA("omp parallel for")
        for(size_t j = 0; j < niter; ++j) {
            const size_t val(1ull << pair.second), hsz(pair.first);
            std::set<std::pair<uint64_t, uint64_t>> lset;
            RNGType gen(mt() + j * (j + 1));
            hll::hlf_t hlf(16, gen(), hsz - 4);
            for(size_t i(0); i < val; ++i) hlf.addh(gen());
            if(std::abs(hlf.report() - val) > (1.03896 / std::sqrt(1ull << pair.first) * val)) {
                lset.emplace(pair);
                //fprintf(stderr, "Warning: Above expected variance for %u, %u.\n", (unsigned)pair.first, (unsigned)pair.second);
            }
            auto diff = hlf.report() - val; // For maximally mixed metaphors
            diffs[j] = diff;
            {
                OMP_PRAGMA("omp critical")
                diffsum += diff;
                ++numlessmore[diff > 0];
                absdiffsum += std::abs(diff);
                sumlessmore[diff > 0] += absdiffsum;
                diff *= diff;
                localnp += (std::abs(hlf.report() - val) <= (1.03896 / std::sqrt(1ull << pair.first) * val));
                gset.insert(lset.begin(), lset.end());
            }
        }
        // auto sqvar = std::sqrt(std::accumulate(std::begin(diffs), std::end(diffs), 0., [div](auto x, auto y) {return x + y / div * y;}));
        {
            //std::fprintf(stdout, "HLF\t%u\t%u\t%lf\t%lf\t%lf\t%zu\t%zu\t%i\t%lf\t%lf\n", (unsigned)pair.first, (unsigned)pair.second, diffsum / div, sqvar, absdiffsum / div, numlessmore[0], numlessmore[1], (int)niter - localnp, sumlessmore[0] / div, sumlessmore[1] / div);
            diffsumsum += diffsum;
        }
        numpass += localnp;
    }
    REQUIRE(numpass >= (niter * pairs.size()));
}

TEST_CASE("phll") {
    spvec_t vec;
    std::vector<size_t> nps {10, 11, 14, 18};
    for(const auto np: nps) {
        const double exact(count_cardinality<score::Lex>(paths, 31, 31, vec, true, nullptr, 2));
        const double inexact(estimate_cardinality<score::Lex>(paths, 31, 31, vec, true, nullptr, 2, np));
        hll::hll_t hll(np);
        {
            Encoder<> enc(Spacer(31, 31, vec), true);
            for(auto &p: paths)
                enc.for_each([&hll](auto x) {hll.addh(x);}, p.data());
            std::fprintf(stderr, "The hll_add method returned an estimated quantity of %lf, compared to estimated %lf manually\n", hll.report(), inexact);
            ssize_t hll_rep(hll.report());
            REQUIRE(hll_rep == size_t(inexact));
        }
        fprintf(stderr, "For np %zu, we have %lf for expected and %lf for measured as correct as we expect to be, theoretically for %lf and %lf.\n",
                np, (1.03896 / std::sqrt(1ull << np)), (std::abs((double)exact - inexact) / inexact), exact, inexact);
        REQUIRE(std::abs(exact - inexact) < hll.est_err() * 2);
    }
}

TEST_CASE("Invert Wang") {
    u64 r;
    std::mt19937_64 mt;
    for(size_t i{0}; i < 10000uL; ++i) {
        REQUIRE(irving_inv_hash(wang_hash(i)) == i);
        r = mt();
        REQUIRE(irving_inv_hash(wang_hash(r)) == r);
        REQUIRE(irving_inv_hash(wang_hash(r)) == r);
    }
}
