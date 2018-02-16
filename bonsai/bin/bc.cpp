#include "util.h"

int main() {
    int c;
    u64 t(0), o(0);
    while((c = fgetc_unlocked(stdin)) != EOF) {
        o += __builtin_popcount(c);
        t += 8;
    }
    std::fprintf(stderr, "1: %zu,0: %zu\n", o, t - o);
}
