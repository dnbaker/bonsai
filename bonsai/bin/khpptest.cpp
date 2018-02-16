#include "getopt.h"
#include "omp.h"
#include <set>
#include <iostream>
#include <fstream>
#include <cassert>
#include <atomic>
#include "khpp.h"

#if 0

int main(int argc, char *argv[]) {
    std::cerr << "Hash map fails all tests. It does not work. Do not use. I will try to entirely rewrite it later.\n";
    kh::khpp<int, int> hello(1 << 20);
    std::fprintf(stderr, "Amount of memory used: %zu\n", hello.estimate_memory());
    return 1;
}
#else
int main() {}
#endif
