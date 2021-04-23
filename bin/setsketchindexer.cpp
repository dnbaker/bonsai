#include "bonsai/encoder.h"
#include "bonsai/util.h"
#include "kseq_declare.h"
#include <queue>
#include "sketch/setsketch.h"
//#include "hll/flat_hash_map/flat_hash_map.hpp"
#include <getopt.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "bonsai/ssi.h"

using namespace bns;

#ifndef CSETFT
#define CSETFT double
#endif

auto gettime() {return std::chrono::high_resolution_clock::now();}

void usage() {
    std::fprintf(stderr, "Usage: setsketchindex <int> <path1> ... <pathN>\nExample: setsketchindex 4096 set1.4096.ss set2.4096.ss\n");
}

struct pq_t: public std::priority_queue<std::pair<double, uint32_t>, std::vector<std::pair<double, uint32_t>>, std::greater<std::pair<double, uint32_t>>> {
    std::vector<std::pair<double, uint32_t>>  &getc() {return this->c;}
    const std::vector<std::pair<double, uint32_t>> &getc() const {return this->c;}
};

int main(int argc, char **argv) {
    int c;
    bool dense_index = false;
    unsigned int k = 5;
    for(;(c = getopt(argc, argv, "k:Dh?")) >= 0;)  {
        switch(c) {
                   case '?': case 'h': usage(); return 1;
                   case 'D': dense_index = true; break;
                   case 'k': k = std::atoi(optarg); break;
        }
    }
    auto t = gettime();
    if(argc == optind || argc + 1 == optind) {usage(); return 1;}
    const size_t m = std::strtoull(argv[optind], nullptr, 10);
    for(auto it = argv + optind + 1; *it; ++it) {
        if(!isfile(*it)) throw std::runtime_error(std::string("File ") + *it + " does not exist.\n");
    }
    bns::SetSketchIndex<CSETFT> ssi(m, dense_index);
    std::vector<CSetSketch<CSETFT>> sketches;
    size_t nfiles = argc - optind - 1;
    sketches.reserve(nfiles);
    for(auto it = argv + optind + 1;*it; ++it) {
        CSetSketch<CSETFT> cs(*it);
        ssi.update(cs);
        sketches.emplace_back(std::move(cs));
    }
    auto t2 = gettime();
    std::fprintf(stderr, "Construction in %gms\n", std::chrono::duration<double, std::milli>(t2 - t).count());
    std::vector<std::pair<std::vector<uint32_t>, std::vector<uint32_t>>> lsh_results(ssi.size());
    const auto ns = ssi.size();
    std::mutex mut;
    t = gettime();
    OMP_PFOR
    for(size_t i = 0; i < ns; ++i) {
        lsh_results[i] = ssi.query_candidates(sketches[i], (k * 3));
        pq_t pq;
        auto &ids = lsh_results[i].first;
        for(size_t L = 0; L < ids.size(); ++L) {
            auto oid = ids[L];
            if(oid == i) continue;
            const double jacc = sketches[i].jaccard_index(sketches[oid]);
            if(pq.size() < k) {
                pq.push({jacc, oid});
            } else if(jacc > pq.top().first) {
                pq.pop();
                pq.push({jacc, oid});
            }
        }
        auto &c = pq.getc();
        std::sort(c.begin(), c.end(), std::greater<void>());
        {
            std::lock_guard<std::mutex> lock(mut);
            std::fprintf(stdout, "%zu", i);
            for(size_t j = 0; j < k; ++j) {
                std::fprintf(stdout, "\t%u:%0.6g", c[j].second, c[j].first);
            }
            std::fputc('\n', stdout);
        }
    }
    t2 = gettime();
    std::fprintf(stderr, "KNN detection in %gms\n", std::chrono::duration<double, std::milli>(t2 - t).count());
    return 0;
}
