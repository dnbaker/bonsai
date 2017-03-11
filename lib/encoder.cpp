#include <encoder.h>

namespace emp {

template std::size_t estimate_cardinality<lex_score>(const std::vector<std::string> &paths,
                                                     unsigned k, uint16_t w, spvec_t spaces,
                                                     void *data, int num_threads, std::size_t np);
template std::size_t estimate_cardinality<hash_score>(const std::vector<std::string> &paths,
                                                     unsigned k, uint16_t w, spvec_t spaces,
                                                     void *data, int num_threads, std::size_t np);

}
