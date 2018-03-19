#define ENABLE_HLL_DEVELOP 1
#include <array>
#include "hll/hll.h"
#include "omp.h"
#include "util.h"
#include <thread>
using namespace emp;
using namespace hll;


std::string calculate_errors(size_t ss, size_t nfiltl2, size_t niter, size_t nelem) {
    uint64_t val;
    std::mt19937_64 mt((ss + 1) * (nfiltl2 + 3) * (niter + 7) * (nelem + 63));
    std::array<double,8> mdiffs {0.,0.,0.,0.,0.,0.};
    if((int)ss - (int)nfiltl2 < 6) {
        std::fprintf(stderr, "Can't do this many.\n");
        return std::string();
    }
    hlf_t hlf(1 << nfiltl2, 137, ss - nfiltl2);
    hll_t hll(ss);
    double frac = 0., fracborrow = 0.;
    for(size_t i(0); i < niter; ++i) {
        for(auto c(nelem); c--;) {
            val = mt();
            hlf.add(val), hll.addh(val);
        }
        mdiffs[0] += std::abs(nelem - hlf.report());
        mdiffs[1] += std::abs(nelem - hll.report());
        mdiffs[2] += std::abs(nelem - hlf.med_report());
        mdiffs[3] += std::abs(nelem - hlf.chunk_report());
        mdiffs[4] += nelem - hlf.report();
        mdiffs[5] += nelem - hll.report();
        mdiffs[6] += nelem - hlf.med_report();
        mdiffs[7] += nelem - hlf.chunk_report();
        frac += nelem / hll.report();
        fracborrow += nelem / hlf.chunk_report();
        hlf.clear();
        hll.clear();
    }
    frac /= niter;
    fracborrow /= niter;
    for(auto &md: mdiffs) md *= 1./niter;
    auto tmp = ks::sprintf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%zu\t%zu\t%zu\n", mdiffs[0], mdiffs[1], mdiffs[2], mdiffs[3],
                           mdiffs[4], mdiffs[5], mdiffs[6], mdiffs[7], frac, fracborrow, ss, size_t(1ull << nfiltl2), nelem);
    return std::string(tmp.begin(), tmp.end());
}

struct comb_t {
    size_t size; size_t nelem; size_t nfiltl2;
};

int main(int argc, char *argv[]) {
    unsigned nthreads = argc > 1 ? std::strtoul(argv[1], nullptr, 10): 8;
    size_t niter      = argc > 2 ? std::strtoull(argv[2], nullptr, 10): 1000;
    std::vector<size_t> sizes    {12, 14, 16, 18, 20};
    std::vector<size_t> nelems   {1 << 18, 1 << 20, 1 << 14, 1 << 12, 1 << 24};
    std::vector<size_t> nfiltl2s {1, 2, 4, 6};
    std::fprintf(stdout, "#Mean error hlf\tMean error hll\tMean error hlf median\t""Mean error strength borrowing\t"
                          "Mean diffs hlf\tMean diffs hll\tMean diffs hlf med\tMean diffs strength borrowing\t"
                          "Mean fraction off (hll)\tmean frac off (hlf borrow)\tsketch size l2\tNumber of subfilters\tnelem\n");
    std::fflush(stdout);
    omp_set_num_threads(nthreads);
    std::string ret;
    ret.reserve(30 * sizes.size() * nelems.size() * nfiltl2s.size());
    std::vector<comb_t> combs;
    for(auto size: sizes) {
        for(auto nelem: nelems) {
            for(auto nfiltl2: nfiltl2s) {
                combs.emplace_back(comb_t{size, nelem, nfiltl2});
            }
        }
    }
    #pragma omp parallel for
    for(unsigned i = 0; i < combs.size(); ++i) {
       auto str = calculate_errors(combs[i].size, combs[i].nfiltl2, niter, combs[i].nelem);
       {
            #pragma omp critical
            ret += str;
       }
    }
    ::write(fileno(stdout), ret.data(), ret.size());
}
