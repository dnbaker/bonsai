#include "bonsai/include/kmeridx.h"

int main(int argc, char *argv[]) {
    if(argc != 3) RUNTIME_ERROR("Usage: <executable> path.fa out.index.gz\n This builds a kmer index from path.fa and writes it to out.index.gz\nIn testing, only k = 31 is supported");
    bns::KmerIdx<uint64_t, uint64_t> idx(31, argv[1]);
    idx.write(argv[2]);
}
