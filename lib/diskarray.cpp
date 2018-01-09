#include "diskarray.h"
namespace ba {
void allocate_file(std::FILE *fp, size_t nbytes) {
    std::fseek(fp, 0, SEEK_SET); // Ensure we're at the beginning.
    static const size_t bufsz(1 << 18);
    std::vector<uint8_t> data(bufsz);
    const int fn(fileno(fp));
    while(nbytes >= bufsz) {
        ::write(fn, data.data(), bufsz);
        nbytes -= bufsz;
    }
    ::write(fn, data.data(), nbytes);
}
}
