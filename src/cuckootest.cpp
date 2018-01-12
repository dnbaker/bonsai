#include "libcuckoo/libcuckoo/cuckoohash_map.hh"
#include "lib/util.h"
#include "lib/encoder.h"
#include "kspp/ks.h"
#include <fstream>
#include "klib/kseq.h"
#include <deque>

#ifndef ITYPE
#define ITYPE u64
#endif

template<typename IntType>
struct identhash {
    __attribute__((always_inline)) std::size_t operator()(IntType val) const {
        return val / sizeof(IntType);
    }
};

void do_work(cuckoohash_map<ITYPE, u64> &ccounter, u64 start, u64 nelem) {
    const auto lambda = [](u64 &num) { ++num; };
    for(ITYPE i(start); i < ITYPE(start + nelem); ++i) {
        ccounter.upsert(i & 255, lambda, 1);
    }
}

using namespace emp;

template<typename TableType>
void update_kmer_table(TableType &ccounter, khash_t(p) *tax, const char *fname, const Spacer &sp, const tax_t taxid) {
    gzFile fp(gzopen(fname, "rb"));
    kseq_t *ks(kseq_init(fp));
    Encoder enc(sp);
    uint64_t kmer;
    auto func = [&](u32 &ctax) {
        ctax = ctax == 0 ? taxid: lca(tax, taxid, ctax);
    };
    while(kseq_read(ks) >= 0) {
        enc.assign(ks);
        while(likely(enc.has_next_kmer())) {
            if((kmer = enc.next_minimizer()) != BF)
                ccounter.upsert(kmer, func, 0);
        }
        if(ccounter.bucket_count() % 256 == 0) std::cerr << "BUCKET COUNT " << ccounter.bucket_count();
    }
    gzclose(fp);
    kseq_destroy(ks);
}

int main(int argc, char *argv[]) {
    if(argc < 3) {
        usage:
        std::fprintf(stderr, "Usage: %s <output_path> <tax_path> <name_table>\n", argv[0]);
        std::exit(1);
    }
    unsigned nthreads(8);
    FILE *ofp = stdout;
    cuckoohash_map<ITYPE, u32, identhash<ITYPE>> ccounter;
    std::deque<std::thread> threads;
    std::string spacing, paths_file;
    int co, k(31), wsz(-1);
    while((co = getopt(argc, argv, "F:w:k:s:c:p:o:h?")) >= 0) {
        switch(co) {
            case 'h': case '?': goto usage;
            case 'p': nthreads = std::atoi(optarg); break;
            case 'F': paths_file = optarg; break;
            case 'o': ofp = std::fopen(optarg, "w"); break;
            case 's': spacing = optarg; break;
            case 'k': k = atoi(optarg); break;
            case 'w': wsz = atoi(optarg); break;
        }
    }
    std::vector<std::string> inpaths;
    if(paths_file.size()) {
        std::ifstream is(paths_file.data());
        std::string line;
        while(std::getline(is, line)) {
            if(line.back() == '\n') line.pop_back();
            if(line.size()) inpaths.emplace_back(std::move(line));
        }
    } else {
        inpaths = std::vector<std::string>(argv + optind + 3, argv + argc);
    }
    if(wsz < 0) wsz = k;
    khash_t(p) *taxmap(build_parent_map(argv[optind + 1]));
    khash_t(name) *name_hash(build_name_hash(argv[optind + 2]));
    spvec_t sv(spacing.size() ? parse_spacing(spacing.data(), k): spvec_t(k - 1, 0));
    Spacer sp(k, wsz, sv);
    for(const auto &path: inpaths) {
        threads.emplace_back(update_kmer_table<decltype(ccounter)>, std::ref(ccounter), taxmap, path.data(), std::cref(sp), get_taxid(path.data(), name_hash));
        if(threads.size() >= nthreads) {
            threads.front().join();
            threads.pop_front();
        }
    }
    for(auto &thread: threads) thread.join();
    auto lt(ccounter.lock_table());
    ks::string ks;
    const int fn(fileno(ofp));
    for(auto &[k, v]: lt) {
        ks.sprintf("key: %s. val: %zu\n", sp.to_string(k).data(), v);
        if(ks.size() & (1 << 16)) ks.write(fn), ks.clear();
    }
    ks.write(ofp), ks.clear();
    std::ofstream out(argv[optind]);
    out << lt;
    if(ofp != stdout) fclose(ofp);
    lt.unlock();
    khash_destroy(taxmap);
    destroy_name_hash(name_hash);
}
