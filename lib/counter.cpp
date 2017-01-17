#include "counter.h"

namespace count {

template class Counter<uint64_t>;

template<typename T>
size_t unique(T &vec, Counter<uint64_t> &c) {
    c.fill(vec);
    return c.size();
}

size_t unique(std::vector<uint64_t> el, Counter<uint64_t> &);

}
