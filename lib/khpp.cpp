#include "lib/khpp.h"
namespace emp {
namespace kh {

template struct khpp_t<std::uint64_t, std::uint64_t, emp::wang_hash_struct>;
using kh64_t = khpp_t<std::uint64_t, std::uint64_t, emp::wang_hash_struct>;
template struct khpp_t<std::uint64_t, std::uint64_t, emp::idt_struct<std::uint64_t>>;

void f() {
    khpp_t<std::uint64_t, std::uint64_t, emp::idt_struct<std::uint64_t>> map;
}

}
}
