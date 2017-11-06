#include <encoder.h>

namespace emp {

template size_t estimate_cardinality<lex_score>(const std::vector<std::string> &paths,
                                                     unsigned k, uint16_t w, spvec_t spaces,
                                                     void *data, int num_threads, size_t np);
template size_t estimate_cardinality<hash_score>(const std::vector<std::string> &paths,
                                                     unsigned k, uint16_t w, spvec_t spaces,
                                                     void *data, int num_threads, size_t np);

}
