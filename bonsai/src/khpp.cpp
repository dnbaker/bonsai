#if 0
#include "khpp.h"
namespace bns {
namespace kh {

template struct khpp_t<u64, u64, bns::wang_hash_struct>;
using kh64_t = khpp_t<u64, u64, bns::wang_hash_struct>;
template struct khpp_t<u64, u64, bns::idt_struct<u64>>;

void f() {
    khpp_t<u64, u64, bns::idt_struct<u64>> map;
}

}
}
#endif
