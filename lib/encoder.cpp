#include "encoder.h"
#include "htslib/kseq.h"
using namespace kpg;

int is_lt(uint64_t i, uint64_t j, void *data) {return i < j;}

int main(void) {
    spvec_t s;
    for(;;) {
        s.push_back(0);
        if(s.size() == 30) break;
        s.push_back(2);
        if(s.size() == 30) break;
        s.push_back(1);
        if(s.size() == 30) break;
    }
    Spacer space(31, 100, &s);
    fprintf(stderr, "13337 in kmer form: %s.\n", space.to_string(13337).data());
    kstring_t ks{0};
    kputs("ACATNGNTNTNATNAAATACACCCCCCCCCCCCTTTTTTTGGGGGGGGGGGGGGGGGN", &ks);
    Encoder<is_lt> enc(ks.s, ks.l, space);
}
