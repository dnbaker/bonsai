#include <unordered_map>
#include <map>
#include <vector>
#include <algorithm>
#include <zlib.h>
#include <stdint.h>
#include <cassert>
#include <string>
#include <fstream>
#include <iostream>
#include "klib/kseq.h"
#include "klib/kthread.h"

KSEQ_INIT(gzFile, gzread)

struct kth {
    const std::string &path_;
    std::map<size_t, size_t> &map_;
    kth(const std::string &path, std::map<size_t, size_t> &map): path_(path), map_(map) {}
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
    size_t nthreads(16);

    gzFile ofp(gzdopen(STDOUT_FILENO, "wT"));

    std::unordered_map<size_t, size_t> data;
    std::vector<std::map<size_t, size_t>> maps;
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
                    if(tmp[tmp.size() - 1] == '\n') tmp.pop_back();
                    paths.push_back(std::move(tmp));
                }
            }
        }
        paths.emplace_back(*p);
    }

    for(size_t i(0); i < paths.size(); ++i) maps.emplace_back();
    std::vector<kth> helpers;
    helpers.reserve(paths.size());
    for(size_t i(0); i < paths.size(); ++i) helpers.emplace_back(paths[i], maps[i]);
    fprintf(stderr, "Processing %zu paths.\n", paths.size());
    kt_for(nthreads, &kfunc, (void *)helpers.data(), paths.size());
    for(auto &kf: helpers)
        for(auto &pair: kf.map_)
            data[pair.first] += pair.second;
    std::vector<size_t> lens;
    lens.reserve(data.size());
    for(auto &pair: data) lens.emplace_back(pair.first);
    std::sort(lens.begin(), lens.end());
    for(auto len: lens)
        gzprintf(ofp, "#Len: %zu\tCount: %zu\n", len, data[len]);
    gzclose(ofp);
}
