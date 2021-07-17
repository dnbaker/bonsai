#ifndef BNS_ALPHBET_H__
#define BNS_ALPHBET_H__
#include <array>
#include <map>
#include <cstdint>
#include <string>

namespace bns {
namespace alph {
using std::size_t;

/*
 *
 * This 
 */

struct Alphabet {
    const char *name;
    const char *setstr;
private:
    size_t nc;
public:
    bool padding;
    size_t nchars() const {return nc + 1;} // One for padding
    using LUType = std::array<int8_t, 256>;
    LUType lut;
    static constexpr LUType make_lut(const char *s, const size_t nc, bool padding=false) {
        LUType arr{-1};
        int id = padding;
        for(size_t ci = 0, i = 0; i < nc; ++i, ++id) {
            while(s[ci] && s[ci] != ',') {
                const auto v = s[ci++];
                arr[v | 32] = arr[v & static_cast<uint8_t>(0xdf)] = id; // lower-case and upper-case
            }
            if(s[ci] == ',') ++ci;
        }
        // Handle Pyrrolysine (P) by mapping it to Lysine (K) if unhandled
        if(arr['P'] == 0) arr['P'] = arr['p'] = arr['K'];
        // Handle SelenoCysteine (U) by mapping it to Cysteine if unhandled
        if(arr['U'] == 0) arr['U'] = arr['u'] = arr['C'];
        return arr;
    }
    constexpr uint8_t translate(char x) const {return lut[static_cast<uint8_t>(x)];}
    static constexpr size_t ncommas(const char *s) {
        size_t ret = 0, i = 0;
        while(s[i]) ret += (s[i] == ','), ++i;
        return ret;
    }
    constexpr Alphabet(const Alphabet &) = default;
    constexpr Alphabet(Alphabet &&) = default;
    constexpr Alphabet(const char *name, const char *s, bool padding=false): name(name), setstr(s), nc(ncommas(s)), padding(padding), lut(make_lut(s, nc, padding)) {
    }
    static constexpr LUType emptylut(bool padding) {
        LUType tlut{-1};
        for(int i = 0; i < 256; ++i) tlut[i] = i + padding;
        return tlut;
    }
    constexpr Alphabet(bool padding=false): name("Bytes"), setstr(""), nc(255 + padding), padding(padding), lut(emptylut(padding)) {
    }
};

// Protein Alphabets
static constexpr const Alphabet BYTES;
static constexpr const Alphabet AMINO20("Standard20", "A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y");

static constexpr const Alphabet SEB14("SE-B(14)", "A,C,D,EQ,FY,G,H,IV,KR,LM,N,P,ST,W");

static constexpr const Alphabet SEB10("SE-B(10)", "AST,C,DN,EQ,FY,G,HW,ILMV,KR,P");
static constexpr const Alphabet SEV10("SE-V(10)", "AST,C,DEN,FY,G,H,ILMV,KQR,P,W");
static constexpr const Alphabet SOLISD("Solis-D", "AM,C,DNS,EKQR,F,GP,HT,IV,LY,W");
static constexpr const Alphabet SOLISG("Solis-G", "AEFIKLMQRVW,C,D,G,H,N,P,S,T,Y");
static constexpr const Alphabet MURPHY("Murphy", "A,C,DENQ,FWY,G,H,ILMV,KR,P,ST");
static constexpr const Alphabet LIA10("Li-A(10)", "AC,DE,FWY,G,HN,IV,KQR,LM,P,ST");
static constexpr const Alphabet LIB10("Li-B(10)", "AST,C,DEQ,FWY,G,HN,IV,KR,LM,P");

static constexpr const Alphabet SEB8("SE-B(8)","AST,C,DHN,EKQR,FWY,G,ILMV,P");
static constexpr const Alphabet SEB6("SE-B(6)","AST,CP,DHNEKQR,FWY,G,ILMV");

static constexpr const Alphabet DAYHOFF("Dayhoff","AGPST,C,DENQ,FWY,HKR,ILMV");

// DNA alphabets

static constexpr const Alphabet DNA4("DNA4", "A,C,G,T");
static constexpr const Alphabet DNA5("DNA5", "A,C,G,T,NMRWSYKVHDB");

static constexpr const Alphabet DNA2KETAMINE("DNA2", "ACM,KGT"); // Amino/Ketones
static constexpr const Alphabet DNA2PYRPUR("DNA2", "AGR,YCT"); // Purines/Pyrimidines


// Source: Reference
// Edgar, RC (2004) Local homology recognition and distance measures in linear time using compressed amino acid alphabets, NAR 32(1), 380-385. doi: 10.1093/nar/gkh180
static const std::map<std::string, const Alphabet *> CAMAP {
    {"AMINO20", &AMINO20},
    {"AMINO", &AMINO20},
    {"PROTEIN", &AMINO20},
    {"SEB8", &SEB8},
    {"SEB10", &SEB10},
    {"SEB14", &SEB14},
    {"SEV10", &SEV10},
    {"MURPHY", &MURPHY},
    {"LIA10", &LIA10},
    {"LIB10", &LIB10},
    {"SEB6", &SEB6},
    {"DAYHOFF", &DAYHOFF},

    {"KETO", &DNA2KETAMINE},
    {"PURPYR", &DNA2PYRPUR},
    {"DNA4", &DNA4},
    {"DNA", &DNA4},
    {"DNA5", &DNA5}
};

} // alph

} // namespace bns


#endif /* BNS_ALPHBET_H__ */
