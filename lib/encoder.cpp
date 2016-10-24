#include "encoder.h"
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
    const Spacer space(31, 100, &s);
    const std::vector<std::string> paths{"/home-1/dbaker49@jhu.edu/work/references/e_coli.fa",
                                         "/home-1/dbaker49@jhu.edu/work/references/lambda_virus.fa"};
    size_t cardinality = estimate_cardinality<is_lt, 18>(paths, 31, 55, nullptr, 0, 2);
    fprintf(stderr, "Number of elements, approximately: %zu.\n", cardinality);
}
