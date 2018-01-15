#pragma once
#include "lib/util.h"
#include "lib/bits.h"
#include <climits>
#include <sys/stat.h>
#include <sys/mman.h>

namespace ba {
using emp::khash_t(64);
using emp::tax_t;

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
        if((mm_ = (char *)mmap(nullptr, sb.st_size, PROT_READ | PROT_WRITE, MAP_SHARED, fileno(fp_), 0)) == MAP_FAILED) throw 1;
        //madvise(mm_, memsz(), MADV_RANDOM);
        if(madvise(mm_, memsz(), MADV_DONTNEED)) throw 1;
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
    void set1_ts(size_t index) {
        for(const char val(1u << (index & 0x7u));
            (mm_[index>>3] & val) == 0;
            __sync_bool_compare_and_swap(mm_ + (index>>3), mm_[index>>3], mm_[index>>3] | val));
    }
    void set0_ts(size_t index) {
        for(const char val(~(1u << (index & 0x7u)));
            mm_[index>>3] & ~val;
            __sync_bool_compare_and_swap(mm_ + (index>>3), mm_[index>>3], mm_[index>>3] & val));
    }
    void set1(size_t index) {
        mm_[index>>3] |= (1u << (index & 0x7u));
    }
    void set0(size_t index) {mm_[index>>3] &= ~(1u << (index & 0x7u));}
    void set1(size_t row, size_t column) {set1(row * nc_ + column);}
    void set0(size_t row, size_t column) {set0(row * nc_ + column);}
    bool operator[](size_t pos) const {return mm_[pos>>3] & (1u << (pos & 0x7u));}
    bool operator()(size_t row, size_t column) const {return operator[](row * nc_ + column);}
    size_t popcount() const {
        size_t ret, i;
        for(i = ret = 0; i < memsz(); ret += emp::popcnt::popcount((unsigned)(uint8_t)mm_[i++]));
        return ret;
    }
    void fill_row(std::vector<uint8_t> &data, size_t row) const {
        if(data.size() != nc_) data.resize(nc_);
        std::copy(data.begin(), data.end(), mm_ + row * nc_);
    }
    std::vector<uint8_t> get_row(size_t row) const {
        std::vector<uint8_t> ret;
        fill_row(ret, row);
        return ret;
    }
};

class MMapTaxonomyBitmap: public DiskBitArray {
    MMapTaxonomyBitmap(size_t nkmers, size_t ntax): DiskBitArray(nkmers, ntax) {}
    void set_kmer(const khash_t(64) *map, u64 kmer, tax_t tax) {
        khint_t ki(kh_get(64, map, kmer));
        if(ki == kh_end(map)) throw 1;
        set1(kh_val(map, ki), tax);
    }
    bool contains_kmer(const khash_t(64) *map, u64 kmer, tax_t tax) const {
        khint_t ki(kh_get(64, map, kmer));
        if(ki == kh_end(map)) throw 1;
        return operator()(kh_val(map, ki), tax);
    }
    void fill_kmer_bitmap(std::vector<uint8_t> &data, const khash_t(64) *map, u64 kmer) {
        if(data.size() != nc_) data.resize(nc_);
    }
    std::vector<uint8_t> kmer_bitmap(const khash_t(64) *map, u64 kmer) {
        std::vector<uint8_t> ret;
        fill_kmer_bitmap(ret, map, kmer);
        return ret;
    }
};

} // namespace ba
