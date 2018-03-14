#include "feature_min.h"
#include "setcmp.h"
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

//#define NOTHREADING

void print_results(ks::string &ks, unsigned sketchsize, double sketchval, double exactval, const hll::hll_t &h1, const khash_t(all) *s1, const hll::hll_t &h2, const khash_t(all) *s2, size_t i, size_t j, unsigned k, char *argv[], int optind) {
    ks.sprintf("%s\t%s\t%lf\t%lf\t%lf\t%lf\t%u\t%u\t%lf\t%zu\n", argv[optind + i], argv[optind + j], sketchval, exactval, std::abs(sketchval - exactval), std::abs(sketchval - exactval) / exactval * 100., sketchsize, k, h1.creport(), kh_size(s1), h2.creport(), kh_size(s2));
}

class cmp_accumulonimbus {
    size_t iherrs_;
    double herrs_;
    double herrsqs_;
    size_t iherrsqs_;
    double uerrs_;
    size_t iuerrs_;
    double uerrsqs_;
    size_t iuerrsqs_;
    double jerrs_;
    double jerrsqs_;
    size_t nadded_;
    size_t npadded_;
public:
    cmp_accumulonimbus() {
        std::memset(this, 0, sizeof(*this));
    }
    void add(const hll::hll_t &h, const khash_t(all) *s) {
        ++nadded_;
        double diff = kh_size(s) - h.creport();
        herrs_ += diff;
        iherrs_ += diff;
        herrsqs_ += diff * diff;
        iuerrs_ += diff * diff;
    }
    void add(const hll::hll_t &h1, const khash_t(all) *s1, const hll::hll_t &h2, const khash_t(all) *s2) {
        ++npadded_;
        auto us = union_size(s1, s2);
        auto est_us = union_size(h1, h2);
        auto diff = us - est_us;
        uerrs_ += diff;
        iuerrs_ += diff;
        diff *= diff;
        uerrsqs_ += diff;
        iuerrsqs_ += diff;
        auto ji = jaccard_index(s1, s2);
        auto est_ji = jaccard_index(h1, h2);
        diff = ji - est_ji;
        jerrs_ += diff;
        diff *= diff;
        jerrsqs_ += diff;
    }
    auto report(std::FILE *fp=stdout) {
        auto ret = std::fprintf(fp, "##total hll errors\tMSE hll\ttotal union errors\tMSE hll\ttotal ji errors\tSE ji\nNum sketches\tNum pairs\n");
        ret     += std::fprintf(fp, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%zu\t%zu\n", herrs_, herrsqs_ / nadded_, uerrs_, uerrsqs_ / npadded_, jerrs_, jerrsqs_, nadded_, npadded_);
    }
};

int main(int argc, char *argv[]) {
    int c;
    bool lowmem(false), canon(true), use_ertl(true), same_stream(false);
    unsigned sketchsize(18), k(31);
    FILE *ofp(stdout), *sumfp(stdout);
    cmp_accumulonimbus accum;
    while((c = getopt(argc, argv, "s:p:o:k:n:mEChS")) >= 0) {
        switch(c) {
            case 'E': use_ertl   = false; break;
            case 'C': canon      = false; break;
            case 'n': sketchsize = std::atoi(optarg); LOG_DEBUG("Set sketch size to %u with string '%s'\n", sketchsize, optarg); break;
            case 'k': k          = std::atoi(optarg); break;
            case 'o': ofp        = std::fopen(optarg, "w"); break;
            case 's': sumfp      = std::fopen(optarg, "w"); break;
            case 'S': same_stream = true;
            case 'p': omp_set_num_threads(std::atoi(optarg)); break;
            case 'm': lowmem     = true; break;
            case 'h': case '?': usage(argv[0]); break;
        }
    }
    if(same_stream) sumfp = ofp;
    std::mutex output_lock;
    omp_set_num_threads(1); // Only using one thread currently as multithreading has not been debugged.
    std::vector<char> buf(1 << 16);
    std::setvbuf(ofp, buf.data(), _IOFBF, buf.size());
    const size_t ngenomes(argc - optind);
    LOG_DEBUG("Comparing %u genomes\n", unsigned(ngenomes));
    if(ngenomes == 0) return EXIT_FAILURE;
    spvec_t sv;
    double sketchval, exactval;
    const int fn(fileno(ofp));
    ks::string ks("#Path1\tPath2\tApproximate jaccard index\tExact jaccard index\t"
                  "Absolute difference\t%difference from exact value\tSketch size\tKmer size\n");
    Spacer sp(k);
    if(lowmem) {
        khash_t(all) *s1(kh_init(all)), *s2(kh_init(all));
        kh_resize(all, s1, 1 << 16), kh_resize(all, s2, 1 << 16);
        hll::hll_t h1(sketchsize, use_ertl), h2(sketchsize, use_ertl);
        std::vector<std::string> paths{"ZOMG"};
        for(size_t i(0); i < ngenomes; ++i) {
            assert(kh_size(s1) == 0);
            LOG_DEBUG("Filling hll\n");
            paths[0] = argv[optind + i];
            fill_set_genome<score::Lex>(argv[optind + i], sp, s1, i, nullptr, canon);
            hll_from_khash(h1, s1);
            h1.sum();
            accum.add(h1, s1);
#ifndef NOTHREADING
            #pragma omp parallel for
#endif
            for(size_t j = i + 1; j < ngenomes; ++j) {
                fill_set_genome<score::Lex>(argv[optind + j], sp, s2, j, nullptr, canon);
                hll_from_khash(h2, s2);
                double sketchval = hll::jaccard_index(h1, h2);
                double exactval  = emp::jaccard_index(s1, s2);
                h2.sum();
                accum.add(h1, s1, h2, s2);
                {
                    #pragma omp critical 
                    print_results(ks, sketchsize, sketchval, exactval, h1, s1, h2, s2, i, j, k, argv, optind);
                }
                kh_clear(all, s2);
                h2.clear();
                if(ks.size() >= 1 << 16)
                {
                    #pragma omp critical 
                    ks.write(fn), ks.clear();
                }
            }
            kh_clear(all, s1);
            h1.clear();
        }
        khash_destroy(s1), khash_destroy(s2);
    } else {
        std::vector<khash_t(all)*> sets;
        while(sets.size() < ngenomes) sets.emplace_back(kh_init(all));
        #pragma omp parallel for
        for(unsigned i = 0; i < ngenomes; ++i) {
            fill_set_genome<score::Lex>(argv[optind + i], sp, sets[i], i, nullptr, canon);
        }
        std::vector<hll::hll_t> sketches;
        while(sketches.size() < ngenomes)
            sketches.emplace_back(
                make_hll(std::vector<std::string>{argv[optind + sketches.size()]},
                         k, k, sv, canon, nullptr, 1, sketchsize, nullptr /*kseq */, use_ertl));
        
        for(unsigned i = 0; i < ngenomes; ++i) {
            accum.add(sketches[i], sets[i]);
        }
        #pragma omp parallel for
        for(size_t i = 0; i < ngenomes; ++i) {
            for(size_t j(i + 1); j < ngenomes; ++j) {
                sketchval = hll::jaccard_index(sketches[i], sketches[j]);
                exactval  = emp::jaccard_index(sets[i], sets[j]);
                {
                    #pragma omp critical
                    accum.add(sketches[i], sets[i], sketches[j], sets[j]);
                    print_results(ks, sketchsize, sketchval, exactval, sketches[i], sets[i], sketches[j], sets[j], i, j, k, argv, optind);
                    if(ks.size() >= 1 << 16) ks.write(fn), ks.clear();
                }
            }
        }
    }
    ks.write(fn); ks.free();
    accum.report(sumfp);
    if(ofp != sumfp && sumfp != stdout) std::fclose(sumfp);
    if(ofp != stdout) std::fclose(ofp);
}
