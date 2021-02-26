#include "sketch/setsketch.h"
#include "zlib.h"

void usage(const char *s) {
    std::fprintf(stderr, "css2ss: Converts a Continuous SetSketch (CSetSketch) to a packed SetSketch.\nUsage: %s <opts>\n-a: set a [required].\n-b: set a [required]\n-B: Pack into bytes\n-N: Pack into nibbles\n-H: Pack into halves\n", s);
    std::exit(1);
}

int main(int argc, char **argv) {
    int c;
    double a = -1., b = -1.;
    int type = 'B';
    while((c = getopt(argc, argv, "HBNh?a:b:")) >= 0) {
        switch(c) {
            case 'h': usage(*argv); return 1;
            case 'a': a = std::atof(optarg); break;
            case 'b': b = std::atof(optarg); break;
            case 'B': type = 'B'; break;
            case 'H': type = 'H'; break;
            case 'N': type = 'N'; break;
        }
    }
    if(a < 0 || b < 0) {
        std::fprintf(stderr, "a, b are mandatory.\n");
        usage(*argv);
    }
    gzFile fp = gzdopen(STDIN_FILENO, "rb");
    sketch::CSetSketch<double> cs(fp);
    gzclose(fp);
    if(type == 'B' || type == 'N') {
        size_t q = (type == 'N' ? size_t(14): size_t(254));
        cs.to_setsketch<uint8_t>(a, b, q).write(stdout);
    } else if(type == 'H') {
        cs.to_setsketch<uint16_t>(a, b).write(stdout);
    }
}
