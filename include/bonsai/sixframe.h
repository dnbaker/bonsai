#ifndef BONSAI_SIXFRAME_H__
#define BONSAI_SIXFRAME_H__
// Six-frame translation of nucleotide sequence for translated (protein)
// k-mer databases, as in Kraken2's translated search: frames +1, +2, +3 on
// the forward strand and -1, -2, -3 on the reverse complement, standard
// genetic code (NCBI table 1). Stop codons split a frame into stop-free
// peptides (k-mers never span a stop). A codon with any base outside ACGT/U
// translates to 'X', which the reduced amino-acid alphabets treat as invalid,
// so k-mers do not span it either.
#include <cstddef>
#include <cstdint>
#include <string>

namespace bns {
namespace sixframe {

// Index 16 * b1 + 4 * b2 + b3 with A = 0, C = 1, G = 2, T/U = 3.
static constexpr char STANDARD_CODE[65] =
    "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF";

inline int base_index(char c) {
    switch(c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': case 'U': case 'u': return 3;
        default: return -1;
    }
}
inline char complement(char c) {
    switch(c) {
        case 'A': case 'a': return 'T';
        case 'C': case 'c': return 'G';
        case 'G': case 'g': return 'C';
        case 'T': case 't': case 'U': case 'u': return 'A';
        default: return 'N';
    }
}
inline char translate_codon(const char *c) {
    const int a = base_index(c[0]), b = base_index(c[1]), d = base_index(c[2]);
    if(a < 0 || b < 0 || d < 0) return 'X';
    return STANDARD_CODE[16 * a + 4 * b + d];
}
// Translation of one frame: s[offset], s[offset + 3], ... (no trailing partial codon)
inline void translate_frame(const char *s, size_t n, unsigned offset, std::string &out) {
    out.clear();
    for(size_t i = offset; i + 3 <= n; i += 3) out.push_back(translate_codon(s + i));
}
inline void reverse_complement(const char *s, size_t n, std::string &out) {
    out.resize(n);
    for(size_t i = 0; i < n; ++i) out[i] = complement(s[n - 1 - i]);
}
// f(peptide pointer, length, frame in {+1,+2,+3,-1,-2,-3}) for every
// non-empty stop-free peptide of every frame, frames in that order.
template<typename F>
void for_each_peptide(const char *s, size_t n, F &&f) {
    std::string rc, pep;
    reverse_complement(s, n, rc);
    for(int strand = 0; strand < 2; ++strand) {
        const char *src = strand ? rc.data(): s;
        for(unsigned off = 0; off < 3; ++off) {
            translate_frame(src, n, off, pep);
            size_t start = 0;
            for(size_t i = 0; i <= pep.size(); ++i) {
                if(i == pep.size() || pep[i] == '*') {
                    if(i > start) f(pep.data() + start, i - start, strand ? -int(off + 1): int(off + 1));
                    start = i + 1;
                }
            }
        }
    }
}

} // namespace sixframe
} // namespace bns

#endif /* BONSAI_SIXFRAME_H__ */
