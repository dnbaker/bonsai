#include "util.h"
#include "getopt.h"

using namespace bns;

int main(int argc, char *argv[]) {
    int c;
    if(argc < 2) {
        usage:
        fprintf(stderr, "Usage: %s <in.txt> out.tax\n" , *argv);
        std::exit(EXIT_FAILURE);
    }
    while((c = getopt(argc, argv, "h?")) >= 0)
        switch(c)
            case 'h': case '?': goto usage;
    khash_t(name) *name_hash(build_name_hash(argv[optind]));
    khash_write<khash_t(name)>(name_hash, argv[optind + 1]);
    destroy_name_hash(name_hash);
    return EXIT_SUCCESS;
}
