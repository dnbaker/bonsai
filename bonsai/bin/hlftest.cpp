#define ENABLE_HLL_DEVELOP 1
#include "hll/hll.h"
#include "util.h"
#include "aesctr.h"
using namespace bns;

int main() {
    aes::AesCtr<std::uint64_t, 8> gen(137);
    std::vector<std::pair<uint64_t, uint64_t>> pairs {
        {10, 14},
        {10, 24},
        {10, 16},
        {10, 18},
        {10, 20},
        {10, 22},
        {10, 24},
        {12, 14},
        {12, 16},
        {12, 18},
        {12, 20},
        {12, 22},
        {12, 24},
        {12, 26},
        {12, 20},
        {15, 21},
        {19, 24},
    };
    int numpass = 0;
    size_t niter = 10000;
    double div = 1./ niter;
    for(const auto &pair: pairs) {
        double diffsum = 0., absdiffsum = 0.;
        std::vector<double> diffs;
        size_t numlessmore[] {0, 0};
        int localnp = 0;
        for(size_t j(0); j < niter; ++j) {
            const size_t val(1ull << pair.second), hsz(pair.first);
            hll::hlf_t hlf(16, gen(), hsz - 4);
            for(size_t i(0); i < val; ++i) {
                hlf.addh(gen());
            }
            if(std::abs(hlf.report() - val) > (1.03896 / std::sqrt(1ull << pair.first) * val)) {
                fprintf(stderr, "Warning: Above expected variance for %u, %u.\n", (unsigned)pair.first, (unsigned)pair.second);
            }
            auto diff = hlf.report() - val; // For maximally mixed metaphors
            diffs.push_back(diff);
            diffsum += diff;
            ++numlessmore[diff > 0];
            absdiffsum += std::abs(diff);
            diff *= diff;
            localnp += (std::abs(hlf.report() - val) <= (1.03896 / std::sqrt(1ull << pair.first) * val));
        }
        auto sqvar = std::sqrt(std::accumulate(std::begin(diffs), std::end(diffs), 0., [div](auto x, auto y) {return x + y / div * y;}));
        std::fprintf(stdout, "HLF\t%u\t%u\t%lf\t%lf\t%lf\t%zu\t%zu\t%i\n", (unsigned)pair.first, (unsigned)pair.second, diffsum / div, sqvar, absdiffsum / div, numlessmore[0], numlessmore[1], (int)niter - localnp);
        numpass += localnp;
    }

}
