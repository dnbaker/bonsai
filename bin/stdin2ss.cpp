#include "sketch/setsketch.h"
void usage() {
    std::fprintf(stderr, "Usage: stdin2ss <flags>\n-S: set sketch size (1000)\n-a: set a\n-b: set b\n-H: use EShortSetS instead of EBytesSetS\n-o: direct output to path instead of stdout\n");
    std::exit(1);
}

int main(int argc, char **argv) {
    int c;
    bool use_short = false;
    double a = -1., b = -1.;
    size_t sketchsz = 1000;
    std::string opath;
    for(;(c = getopt(argc, argv, "Hh?S:a:o:b:")) >= 0;) {
        switch(c) {
            case 'a': a = std::atof(optarg); break;
            case 'b': b = std::atof(optarg); break;
            case 'S': sketchsz = std::strtoull(optarg, nullptr, 10); break;
            case 'H': use_short = true; break;
            case 'o': opath = optarg; break;
            case 'h': case '?': usage();
        }
    }
    std::FILE *ifp = argc == optind ? stdin: std::fopen(argv[optind], "rb");
    uint64_t buf[2];
    if(use_short) {
        if(a < 0) a = sketch::setsketch::EShortSetS::DEFAULT_A;
        if(b < 0) b = sketch::setsketch::EShortSetS::DEFAULT_B;
        sketch::setsketch::EShortSetS ss(sketchsz, b, a);
        while(std::fread(buf, sizeof(buf), 1, ifp) == 1) {
            ss.add(buf[0]);
        }
    } else {
        std::FILE *ofp;
        if(opath.size()) ofp = std::fopen(opath.data(), "wb");
        else ofp = stdout;
        std::fprintf(stderr, "now reading!\n");
        if(a < 0) a = sketch::setsketch::EByteSetS::DEFAULT_A;
        if(b < 0) b = sketch::setsketch::EByteSetS::DEFAULT_B;
        sketch::setsketch::EByteSetS ss(sketchsz, b, a);
        size_t i = 0;
        size_t fret;
        while((fret = std::fread(buf, sizeof(buf), 1, ifp)) == 1u) {
            ++i;
            ss.add(buf[0]);
        }
        ss.write(ofp);
        if(ofp != stdout) std::fclose(ofp);
        std::fprintf(stderr, "processed %zu total\n", i);
    }
    if(ifp != stdin) std::fclose(ifp);
}
