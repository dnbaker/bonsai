#include "sketch/setsketch.h"
#include "sketch/hll.h"
#include "sketch/rng.h"

using sketch::CSetSketch;
using sketch::hll_t;

void usage() __attribute__((noreturn)) {
    std::fprintf(stderr, "errexp <flags>\n-l <linspace = .1>. Set the spacing between 0 and 1 for Jaccard's to measure.\n");
}


int main(int argc, char **argv) {
    std::vector<size_t> skz {{6, 8, 10, 12, 14, 16}}; // sketch sizes
    std::vector<size_t> ssz {{8, 12, 16, 20, 24, 28}}; // set sizes
    double linspace = .1;
    for(int c;(c = getopt(argc, argv, "l:h?")) >= 0;) {
        switch(c) {
            default: case '?': case 'h': usage();
            case 'h': linspace = std::atof(optarg); if(linspace > 1. || linspace < 0.) throw std::invalid_argument("linspace must be > 0 and < 1");
        }
    }
    std::transform(ssz.begin(), ssz.end(), ssz.begin(), [](auto x) {return 1ull << x;});
    std::vector<double> jaccards;
    for(double li = linspace; li < 1.; li = std::min(std::max(li + linspace, li), 1.))
        jaccards.push_back(li);

    for(const size_t setsize: ssz) {
        const size_t ns = skz.size() * jaccards.size();
        std::vector<hll_t> hlls; hlls.reserve(skz.size());
        std::vector<CSetSketch<double>> css;
        for(const auto v: skz) {
            hlls.emplace_back(v);
            css.emplace_back(size_t(1) << v);
        }
        std::vector<hll_t> fhlls = hlls;
        std::vector<CSetSketch<double>> fcss = css;
        for(const auto ji: jaccards) {
            for(auto &h: hlls) {
                h.reset();
            }
            for(auto &h: fhlls) h.reset();
            for(auto &c: css) c.reset();
            for(auto &c: fcss) c.reset();
            size_t i, e = std::ceil(ji * setsize);
            for(auto &h: hlls) {
                for(i = 0; i < e; ++i)
                    h.add(i);
                for(; i < setsize; ++i)
                    h.add(i + 0xfffffffull);
            }
            for(auto &h: css) {
                for(i = 0; i < e; ++i)
                    h.add(i);
                for(; i < setsize; ++i)
                    h.add(i + 0xfffffffull);
            }
            for(auto &h: flls) {
                for(size_t i = 0; i < setsize; ++i)
                    h.add(i);
            }
            double exact_j = double(e) / setsize;
            auto hit = hlls.begin(), fhit = fhlls.begin();
            auto chit = css.begin(), fchit = fcss.begin();
            for(const auto sz: skz) {
                double hji = hit->jaccard_index(*fhit);
                double cji
            }
        }
        }
    }
}
