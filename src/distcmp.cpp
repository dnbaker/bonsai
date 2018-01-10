#include "lib/feature_min.h"
#include "lib/setcmp.h"
#include <omp.h>
#include <getopt.h>

using namespace emp;

void usage(const char *arg) {
    std::fprintf(stderr, "Usage: %s <opts> [genome paths]\n"
                         "-k\tkmer size [31]\n"
                         "-n\tSketch size [18]\n"
                         "-o\tOutput path [stdout]\n"
                         "-m\tUse more CPU but less memory [default: more memory, less runtime]\n"
                         "-h\tEmit usage\n",
                 arg);
    std::exit(EXIT_FAILURE);
}

#if !NOTHREADING
std::mutex output_lock;
#endif

#if NOTHREADING
#define EMIT_RESULTS(sketchval, exactval) do { \
    std::fprintf(ofp, "%s\t%s\t%s\t%lf\t%lf\t%lf\t%lf\n", \
                 argv[optind + i], argv[optind + j], BOOLS[i == j], \
                 sketchval, exactval, std::abs(sketchval - exactval), \
                 std::abs(sketchval - exactval) / exactval * 100.); } while(0)
#else
#define EMIT_RESULTS(sketchval, exactval) do { std::lock_guard<std::mutex> lock(output_lock);\
    std::fprintf(ofp, "%s\t%s\t%s\t%lf\t%lf\t%lf\t%lf\n", \
                 argv[optind + i], argv[optind + j], BOOLS[i == j], \
                 sketchval, exactval, std::abs(sketchval - exactval), \
                 std::abs(sketchval - exactval) / exactval * 100.); } while(0)
#endif

int main(int argc, char *argv[]) {
    int c;
    bool lowmem(false);
    unsigned sketchsize(18), k(31);
    FILE *ofp(stdout);
    while((c = getopt(argc, argv, "p:o:k:n:mh")) >= 0) {
        switch(c) {
            case 'n': sketchsize = std::atoi(optarg); break;
            case 'k':          k = std::atoi(optarg); break;
            case 'h': case '?': usage(argv[0]); break;
            case 'o':       ofp = std::fopen(optarg, "w"); break;
            case 'p': omp_set_num_threads(std::atoi(optarg)); break;
            case 'm': lowmem    = true; break;
        }
    }
    omp_set_num_threads(1); // Only using one thread currently as multithreading has not been debugged.
    static const char *BOOLS [2] {
        "False", "True"
    };
    std::vector<char> buf(1 << 16);
    std::setvbuf(ofp, buf.data(), _IOFBF, buf.size());
    const size_t ngenomes(argc - optind);
    LOG_DEBUG("Comparing %u genomes\n", unsigned(ngenomes));
    spvec_t sv;
    double sketchval, exactval;
    std::fprintf(ofp, "#Path1\tPath2\tTrue if path1==path2\tApproximate jaccard index\tExact jaccard index\t"
                      "Absolute difference\t%% difference from exact value\n");
#if 0
    std::mt19937_64 mt(std::time(nullptr));
    std::shuffle(argv + optind, argv + argc, mt);
#endif
    if(lowmem) {
        khash_t(all) *s1(kh_init(all)), *s2(kh_init(all));
        hll::hll_t h1(sketchsize), h2(sketchsize);
        Spacer sp(k);
        for(size_t i(0); i < ngenomes; ++i) {
            assert(kh_size(s1) == 0);
            fill_hll(h1, std::vector<std::string>{argv[optind + i]}, k, k, sv, nullptr, 1, sketchsize);
            fill_set_genome<score::Lex>(argv[optind + i], sp, s1, i, nullptr);
            #pragma omp parallel for
            for(size_t j = i + 1; j < ngenomes; ++j) {
                fill_set_genome<score::Lex>(argv[optind + j], sp, s2, j, nullptr);
#if COPY_HLL
                h2 = make_hll(std::vector<std::string>{argv[optind + j]}, k, k, sv, nullptr, 1, sketchsize);
#else
                fill_hll(h2, std::vector<std::string>{argv[optind + j]}, k, k, sv, nullptr, 1, sketchsize);
#endif
                double sketchval = hll::jaccard_index(h1, h2);
                double exactval  = emp::jaccard_index(s1, s2);
                EMIT_RESULTS(sketchval, exactval);
                kh_clear(all, s2);
#if COPY_HLL
#else
                h2.clear();
#endif
            }
            kh_clear(all, s1);
            h1.clear();
        }
        khash_destroy(s1), khash_destroy(s2);
    } else {
        std::vector<khash_t(all)*> sets;
        while(sets.size() < ngenomes) sets.emplace_back(kh_init(all));
        Spacer sp(k);
        #pragma omp parallel for
        for(unsigned i = 0; i < ngenomes; ++i) {
            fill_set_genome<score::Lex>(argv[optind + i], sp, sets[i], i, nullptr);
        }
        std::vector<hll::hll_t> sketches;
        while(sketches.size() < ngenomes)
            sketches.emplace_back(
                make_hll(std::vector<std::string>{argv[optind + sketches.size()]},
                         k, k, sv, nullptr, 1, sketchsize));
        for(size_t i(0); i < ngenomes; ++i) {
            for(size_t j(i + 1); j < ngenomes; ++j) {
                sketchval = hll::jaccard_index(sketches[i], sketches[j]);
                exactval  = emp::jaccard_index(sets[i], sets[j]);
                EMIT_RESULTS(sketchval, exactval);
            }
        }
    }
    if(ofp != stdout) std::fclose(ofp);
}
