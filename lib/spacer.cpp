#include "spacer.h"
#include "encoder.h"

// Simple test for spacer IO.
using namespace kpg;

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
