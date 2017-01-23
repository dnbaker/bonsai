#include "spacer.h"
#include "encoder.h"

// Simple test for spacer IO.
namespace emp {


std::uint32_t comb_size(const spvec_t &spaces) {
    std::uint32_t ret(spaces.size() + 1); // Since there's 1 fewer entry in spaces
    // We then increment the size of our comb for each space.
    for(auto i: spaces) {
        ret += i;
    }
    return ret;
}

spvec_t parse_spacing(const char *ss, unsigned k) {
    if(!ss) return spvec_t(k - 1, 0); // No spaces

    spvec_t ret;
    while(ss && *ss) {
        int j(atoi(ss));
        ret.emplace_back(j);

        if(strchr(ss, 'x')) {
            ss = strchr(ss, 'x') + 1;
            for(int k(atoi(ss) - 1); k; k--) ret.emplace_back(j);
        }
        ss = strchr(ss, ',') + 1;
    }
    return ret;
}

}
#ifdef __SPVEC_MAIN
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
}
#endif // #ifdef __SPVEC_MAIN
