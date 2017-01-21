#include "cms.h"

using namespace emp;

#ifdef __CMS_MAIN
int main(void) {
    cms_t<uint8_t> test(65536);
    for(std::size_t i(0);i<100000; ++i) test.add(i);
    fprintf(stderr, "Estimate %u.\n", (unsigned)test.query(1337));
}
#endif
