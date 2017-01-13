#include "counter.h"

namespace count {

template class Counter<uint64_t>;

size_t unique(std::vector<uint64_t> &vec, Counter<uint64_t> c) {
    c.fill(vec);
    return c.size();
}

}
