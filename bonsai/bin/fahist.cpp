#include <unordered_map>
#include <map>
#include <vector>
#include <algorithm>
#include <stdint.h>
#include <cassert>
#include <string>
#include <fstream>
#include <iostream>
#include "klib/kseq.h"
#include "klib/kthread.h"
#if ZWRAP_USE_ZSTD
#  include "zstd_zlibwrapper.h"
#else
#  include <zlib.h>
#endif

#ifdef USE_PDQSORT
#include "pdqsort/pdqsort.h"
#ifndef SORT
#define SORT pdqsort
#endif
#define SORT_BRANCHLESS pdqsort_branchless
#else
#ifndef SORT
#define SORT ::std::sort
#endif
#define SORT_BRANCHLESS ::std::sort
#endif


KSEQ_INIT(gzFile, gzread)

struct kth {
    const std::string &path_;
    std::map<std::size_t, std::size_t> &map_;
    kth(const std::string &path, std::map<std::size_t, std::size_t> &map): path_(path), map_(map) {}
};

void kfunc(void *data_, long index, int id) {
    kth *data(((kth *)data_) + index);
    gzFile ifp(gzopen(data->path_.data(), "r"));
    kseq_t *ks(kseq_init(ifp));
    while(kseq_read(ks) >= 0) ++data->map_[ks->seq.l];
    gzclose(ifp);
    kseq_destroy(ks);
}

int main(int argc, char **argv) {
    std::size_t nthreads(16);

    gzFile ofp(gzdopen(STDOUT_FILENO, "wT"));

    std::unordered_map<std::size_t, std::size_t> data;
    std::vector<std::map<std::size_t, std::size_t>> maps;
    std::vector<std::string> paths;
    paths.reserve(argc - 1);
    
    for(char **p(argv + 1); *p; ++p) {
        if(**p == '-') {
            if((*p)[1] == 'f') {
                ++p;
                paths.clear();
                std::string tmp;
                std::ifstream ifs(*p);
                while(std::getline(ifs, tmp)) {
                    std::cout << tmp << '\n';
                    paths.emplace_back(std::move(tmp));
                }
            }
        }
        paths.emplace_back(*p);
    }

    for(std::size_t i(0); i < paths.size(); ++i) maps.emplace_back();
    std::vector<kth> helpers;
    helpers.reserve(paths.size());
    for(std::size_t i(0); i < paths.size(); ++i) helpers.emplace_back(paths[i], maps[i]);
    fprintf(stderr, "Processing %zu paths.\n", paths.size());
    kt_for(nthreads, &kfunc, (void *)helpers.data(), paths.size());
    for(auto &kf: helpers)
        for(auto &pair: kf.map_)
            data[pair.first] += pair.second;
    std::vector<std::size_t> lens;
    lens.reserve(data.size());
    for(const auto &pair: data) lens.emplace_back(pair.first);
    SORT(lens.begin(), lens.end());
    for(const auto len: lens)
        gzprintf(ofp, "#Len: %zu\tCount: %zu\n", len, data[len]);
    gzclose(ofp);
}
