#include <cstdio>
#include <utility>
#include <cassert>
#include <string>
#include <fstream>
#include <algorithm>
#include <vector>
#include <getopt.h>
#include <omp.h>
#include <zlib.h>
#include "include/sketch/common.h"

void usage() {
    std::fprintf(stderr, "Usage: cmpshs -F fnames.txt > out.bin\n");
    std::exit(1);
}

template<typename T>
static std::vector<T> load_gs(std::string path, int res=1<<16) {
    std::vector<T> ret;
    ret.reserve(res);
    gzFile fp = gzopen(path.data(), "rb");
    T hv;
    while(gzread(fp, &hv, sizeof(hv)) == static_cast<typename std::make_signed<T>::type>(hv)) {
        ret.push_back(hv);
    }
    gzclose(fp);
    return ret;
}

// See https://stackoverflow.com/questions/32640327/how-to-compute-the-size-of-an-intersection-of-two-stl-sets-in-c
struct Counter
{
  struct value_type { template<typename T> value_type(const T&) { } };
  void push_back(const value_type&) { ++count; }
  size_t count = 0;
};


template<typename T, typename A1, typename A2>
size_t intersection_size(const std::vector<T, A1> &a, const std::vector<T, A2> &b) {
    Counter c;
    std::set_intersection(a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(c));
    return c.count;
}


static inline size_t ij2ind(size_t i, size_t j, size_t n) {
    return (((i) * (n * 2 - i - 1)) / 2 + j - (i + 1));
}

template<typename T>
void perform_work(const std::vector<std::string> &paths, std::string opath);

int main(int argc, char *argv[]) {
    int c, nt = 1;
    bool use_32bit = false, use_16bit = false;
    std::string inpaths;
    std::string opath;
    while((c = getopt(argc, argv, "h36p:F:o:")) >= 0) {
        switch(c) {
            case 'p': nt = std::atoi(optarg); break;
            case 'F': inpaths = optarg; break;
            case '3': use_32bit = true; break;
            case '6': use_16bit = true; break;
            case 'o': opath = optarg; break;
            default: case 'h': usage();
        }
    }
    if(use_32bit && use_16bit) throw 7;
    nt = std::max(nt, 1);
    omp_set_num_threads(nt);
    std::vector<std::string> paths;
    std::ifstream ifs(inpaths.size() ? inpaths.data() : argc == optind ? "/dev/stdin": (const char *)argv[optind]);
    for(std::string s; std::getline(ifs, s);paths.push_back(s));
    if(paths.empty()) usage();
    auto fp = use_32bit ? perform_work<uint32_t>
                        : use_16bit ? perform_work<uint16_t>
                                    : perform_work<uint64_t>;
    fp(paths, opath);
}

template<typename T>
void perform_work(const std::vector<std::string> &paths, std::string opath) {
    std::vector<std::vector<T>> genome_sets;
    size_t n = paths.size();
    genome_sets.resize(n);
    OMP_PRAGMA("omp parallel for")
    for(size_t i = 0; i < n; ++i) {
        genome_sets[i] = load_gs<T>(paths[i]);
        assert(std::is_sorted(genome_sets[i].begin(), genome_sets[i].end()));
    }
    //std::vector<std::pair<float, uint32_t>> ret(n * (n - 1) / 2);
    std::vector<float> ret(n * (n - 1) / 2);
    for(size_t i = 0; i < n; ++i) {
        const auto &gs1 = genome_sets[i];
        OMP_PRAGMA("omp parallel for")
        for(size_t j = i + 1; j < n; ++j) {
            auto isz = intersection_size(gs1, genome_sets[j]);
            ret[ij2ind(i, j, n)] = gs1.size() + genome_sets[j].size() - isz;
        }
    }
    if(opath.empty()) opath = "/dev/stdout";
    gzFile fp = gzopen(opath.data(), "wb");
    if(!fp) throw 2;
    if(gzwrite(fp, &n, sizeof(n)) != sizeof(n)) throw 3;
    if(gzwrite(fp, ret.data(), ret.size() * sizeof(ret[0])) != ssize_t(ret.size() * sizeof(ret[0]))) throw 4;
    gzclose(fp);
}
