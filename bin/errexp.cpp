#undef NDEBUG
#include "sketch/setsketch.h"
#include "sketch/hll.h"

using sketch::CSetSketch;
using sketch::hll_t;

void usage()  {
    std::fprintf(stderr, "errexp <flags>\n-l <linspace = .1>. Set the spacing between 0 and 1 for Jaccard's to measure. Emits max/min register values to stderr and a table to stdout\n");
}


int main(int argc, char **argv) {
    const std::vector<size_t> skz {{8, 9, 10, 11, 12, 14}}; // sketch sizes
    std::vector<size_t> ssz {{16, 20, 22, 24, 28}}; // set sizes
    double linspace = .1;
    for(int c;(c = getopt(argc, argv, "l:h?")) >= 0;) {
        switch(c) {
            default: case '?': case 'h': usage(); std::exit(1); break;
            case 'l': linspace = std::atof(optarg); if(linspace > 1. || linspace < 0.) throw std::invalid_argument("linspace must be > 0 and < 1");
        }
    }
    std::transform(ssz.begin(), ssz.end(), ssz.begin(), [](auto x) {return 1ull << x;});
    std::vector<double> jaccards;
    for(double li = linspace; li < 1.; li = std::min(std::max(li + linspace, li), 1.))
        jaccards.push_back(li);

    for(const size_t setsize: ssz) {
        std::vector<hll_t> hlls, fhlls; hlls.reserve(skz.size());
        std::vector<CSetSketch<double>> css, fcss;
        for(const auto v: skz) {
            hlls.emplace_back(v);
            fhlls.emplace_back(v);
            css.emplace_back(size_t(1) << v);
            fcss.emplace_back(size_t(1) << v);
        }
        for(const auto &ji: jaccards) {
            if(&ji != jaccards.data()) {
                for(auto &h: hlls) h.reset();
                for(auto &h: fhlls) h.reset();
                for(auto &c: css) c.reset();
                for(auto &c: fcss) c.reset();
            }
            size_t i, e = std::min(size_t(std::ceil(ji * setsize)), setsize);
            for(auto &h: hlls) {
                for(i = 0; i < e; ++i)
                    h.addh(i);
                for(; i < setsize; ++i)
                    h.addh(i + 0xfffffffull);
                h.sum();
            }
            for(auto &h: css) {
                for(i = 0; i < e; ++i)
                    h.addh(i);
                for(; i < setsize; ++i)
                    h.addh(i + 0xfffffffull);
            }
            for(auto &h: fcss) {
                for(size_t i = 0; i < setsize; ++i)
                    h.addh(i);
            }
            for(auto &h: fhlls) {
                for(size_t i = 0; i < setsize; ++i)
                    h.addh(i);
                h.sum();
            }
            double exact_j = double(e) / (2 * setsize - e);
            auto hit = hlls.begin(), fhit = fhlls.begin();
            auto chit = css.begin(), fchit = fcss.begin();
            assert(hlls.size() == skz.size());
            assert(fhlls.size() == skz.size());
            assert(css.size() == skz.size());
            assert(fcss.size() == skz.size());
            for(const auto sz: skz) {
                if(hit == hlls.end()) throw std::runtime_error("hit is at the end");
                if(fhit == fhlls.end()) throw std::runtime_error("fhit is at the end");
                double hji = hit->jaccard_index(*fhit);
                double cji = chit->jaccard_index(*fchit);
                std::fprintf(stdout, "HLL\t%zu\t%zu\t%0.12g\t%0.12g\t%0.12g\n", sz, setsize, hji, exact_j, exact_j - hji);
                std::fprintf(stdout, "CSSDouble\t%zu\t%zu\t%0.12g\t%0.12g\t%0.12g\n", sz, setsize, cji, exact_j, exact_j - cji);
                auto lhc = chit->cardinality(), rhc = fchit->cardinality();
                {
                    auto lhs = chit->to_setsketch<uint16_t>(1.0006, .001);
                    auto rhs = fchit->to_setsketch<uint16_t>(1.0006, .001);
                    std::fprintf(stderr, "lhs max, min %ld, %ld. rhs %ld/%ld\n", lhs.max(), lhs.min(), rhs.max(), rhs.min());
                    auto abmu16 = lhs.alpha_beta_mu(rhs, lhc, rhc);
                    auto a16 = std::get<0>(abmu16);
                    auto b16 = std::get<1>(abmu16);
                    //auto mu16 = std::get<2>(abmu16);
                    double jac16 = std::max(0., double(1. - a16 - b16));
                    std::fprintf(stdout, "CSSu16\t%zu\t%zu\t%0.12g\t%0.12g\t%0.12g\n", sz, setsize, jac16, exact_j, exact_j - jac16);
                }
                {
                    auto lhs = chit->to_setsketch<uint8_t>(1.11, 100);
                    auto rhs = fchit->to_setsketch<uint8_t>(1.11, 100);
                    std::fprintf(stderr, "lhs max, min %ld, %ld. rhs %ld/%ld\n", lhs.max(), lhs.min(), rhs.max(), rhs.min());
                    auto abmu8 = lhs.alpha_beta_mu(rhs, lhc, rhc);
                    auto a8 = std::get<0>(abmu8);
                    auto b8 = std::get<1>(abmu8);
                    //auto mu8 = std::get<2>(abmu8);
                    double jac8 = std::max(0., double(1. - a8 - b8));
                    std::fprintf(stdout, "CSSu8\t%zu\t%zu\t%0.12g\t%0.12g\t%0.12g\n", sz, setsize, jac8, exact_j, exact_j - jac8);
                }
                ++hit, ++chit, ++fhit, ++fchit;
            }
        }
    }
}
