#if 0
#include "khpp.h"
namespace emp {
namespace kh {

template struct khpp_t<u64, u64, emp::wang_hash_struct>;
using kh64_t = khpp_t<u64, u64, emp::wang_hash_struct>;
template struct khpp_t<u64, u64, emp::idt_struct<u64>>;

void f() {
    khpp_t<u64, u64, emp::idt_struct<u64>> map;
}

}
}
#endif
