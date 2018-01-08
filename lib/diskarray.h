#pragma once
#include "lib/util.h"
#include "lib/bits.h"
#include <climits>
#include <sys/stat.h>
#include <sys/mman.h>

namespace ba {

void allocate_file(std::FILE *fp, size_t nbits);

class DiskBitArray {
protected:
    std::FILE *fp_; // file pointer
    char      *mm_; // mmap pointer
    const size_t nr_, nc_; // Number of rows, columns
    std::string fpath_;
public:
    DiskBitArray(size_t nfeat, size_t nclades, std::string path="./__file.mm"):
        fp_(std::fopen(path.data(), "w+")), mm_(nullptr), nr_(nfeat), nc_(nclades), fpath_(std::move(path)) {
        if(fp_ == nullptr) throw 1;
        allocate_file(fp_, memsz());
        struct stat sb;
        fstat(fileno(fp_), &sb);
        std::fprintf(stderr, "Size of file: %zu\n", size_t(sb.st_size));
        if((mm_ = (char *)mmap(nullptr, sb.st_size, PROT_READ | PROT_WRITE, MAP_SHARED, fileno(fp_), 0)) == MAP_FAILED) throw 1;
        std::fprintf(stderr, "mm_: %p\n", (void *)mm_);
    }
    size_t memsz() const {
        size_t ret(nr_ * nc_);
        ret = (ret >> 3) + !!(ret & 0x7u);
        return ret;
    }
    size_t size() const {return nr_ * nc_;}
    ~DiskBitArray() {
        if(fp_) std::fclose(fp_);
        if(mm_) munmap((void *)mm_, memsz());
    }
    void set1(size_t index) {
        mm_[index>>3] |= (1u << (index & 0x7u));
    }
    void set0(size_t index) {
        mm_[index>>3] &= ~(1u << (index & 0x7u));
    }
    void set1(size_t row, size_t column) {set1(row * nc_ + column);}
    void set0(size_t row, size_t column) {set0(row * nc_ + column);}
    bool operator[](size_t pos) const {return mm_[pos>>3] & (1u << (pos & 0x7u));}
    bool operator()(size_t row, size_t column) const {return operator[](row * nc_ + column);}
    size_t popcount() const {
        using namespace emp;
        u64 *ptr((u64 *)mm_);
        size_t ret(0), i(0);
        while(i < memsz() >> 3) ret += popcnt::popcount(ptr[i++]);
        for(i = memsz() & ~0x7uLL; i < memsz(); ret += popcnt::popcount(mm_[i++]));
        return ret;
    }
};

} // namespace ba
