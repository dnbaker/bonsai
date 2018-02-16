#include "diskarray.h"
#include "aesctr.h"

namespace ba {
void allocate_file(std::FILE *fp, size_t nbytes) {
    std::fseek(fp, 0, SEEK_END);
#if !NDEBUG
    auto npos(std::ftell(fp));
    std::fseek(fp, 0, SEEK_SET); // Ensure we're at the beginning.
    LOG_DEBUG("Size of file before allocation: %s\n", std::to_string(npos).data());
#endif
    static const size_t bufsz(1 << 18);
    std::vector<uint8_t> data(bufsz);
    const int fn(fileno(fp));
    aes::AesCtr<u64> gen(nbytes);
    while(nbytes >= bufsz) {
        ::write(fn, data.data(), bufsz);
        nbytes -= bufsz;
        for(size_t i(0); i < bufsz / sizeof(u64); ++i)
            *((uint64_t *)data.data() + i) = gen();
    }
    ::write(fn, data.data(), nbytes);
}
}
