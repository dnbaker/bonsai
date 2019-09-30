#include "hll/include/sketch.h"
#include "bonsai/include/encoder.h"
#include "omp.h"

using namespace bns;

void usage() {
    std::fprintf(stderr, "rolling_multk_sketch <opts> in.fa\n-P: set prefix\n-p: set number of threads\n-k: Add kmer length\n-C: Do not canonicalize\n-r:  set START,END for kmer range (e.g., -r34,38 will use 34, 35, 36, 37]).\n");
    std::fprintf(stderr, "-S: set log 2 sketch size (14)\n");
    std::exit(1);
}

template<typename IT=uint64_t, typename ... Args>
void build_multiset(IT k, std::vector<gzFile> &fps, std::vector<gzFile> &ofps, bool canon=false, Args &&... args) {
    #pragma omp parallel for
    for(size_t i = 0; i < fps.size(); ++i) {
        auto fp = fps[i];
        bns::Encoder<> enc(k, canon);
        khash_t(p) kh{0,0,0,0,0,0,0};

        // Load all kmers
        enc.for_each([p=&kh](IT hashvalue) {
            hashvalue = sketch::WangHash()(hashvalue); // Use WangHash
            khint_t ki;
            if((ki = kh_get(p, p, hashvalue)) == kh_end(p)) {
                int khr;
                ki = kh_put(p, p, hashvalue, &khr);
                kh_val(p, ki) = 1;
            } else ++kh_val(p, ki);
        }, fp);
        gzclose(fp);
        fp = ofps[i]; // Note reuse of fp for outfile

        // Sort kmers with counts by kmer
        std::vector<std::pair<uint64_t, uint32_t>> tmp;
        tmp.reserve(kh_size(&kh));
        for(khint_t ki = 0; ki != kh_end(&kh); ++ki) {
            if(kh_exist(&kh, ki))
                tmp.push_back(std::make_pair(kh_key(&kh, ki), kh_val(&kh, ki)));
        }
        SORT(tmp.begin(), tmp.end(), [](auto x, auto y) {return x.first > y.first;});
        // Reverse-sorted, like default in sketch library
        // Get summary infomation, write to file
        size_t sum = 0;
        double sqsum = 0;
        for(const auto &p: tmp)
            sum += p.second, sqsum += p.second * p.second;
        uint64_t nelem = tmp.size();
        gzwrite(fp, &nelem, sizeof(nelem));
        gzwrite(fp, &sum, sizeof(sum));
        gzwrite(fp, &sqsum, sizeof(sqsum));
        for(const auto &p: tmp)
            gzwrite(fp, &p.first, sizeof(p.first));
        for(const auto &p: tmp)
            gzwrite(fp, &p.second, sizeof(p.second));
        // Clean up
        gzclose(fp);
        std::free(kh.keys);
        std::free(kh.flags);
        std::free(kh.vals);
    }
}


int main(int argc, char *argv[]) {
    int c;
    std::string prefix;
    int canon = true;
    int k = 31;
    while((c = getopt(argc, argv, "ChP:p:k:r:")) >= 0) {
        switch(c) {
            case 'h': usage(); break;
            case 'k': k = std::atoi(optarg); break;
            case 'C': canon = false; break;
            case 'P': prefix = optarg; break;
            case 'p': omp_set_num_threads(std::atoi(optarg)); break;
        }
    }
    if(optind == argc) usage();
    std::vector<gzFile> fps;
    std::vector<gzFile> ofps;
    if(prefix.empty()) prefix = argv[1];
    for(int i = optind; i < argc; ++i) {
        fps.push_back(gzopen(argv[i], "rb"));
        ofps.push_back(gzopen((std::string(argv[i]) + ".multiset").data(), "wb"));
    }
    build_multiset(k, fps, ofps, canon);
    for(int i = optind; i < argc; ++i) {
        mh::FinalCRMinHash<uint64_t, std::greater<>, uint32_t> ms((std::string(argv[i]) + ".multiset").data());
    }
}
