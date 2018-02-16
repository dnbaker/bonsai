#include <encoder.h>

namespace emp {

template size_t estimate_cardinality<score::Lex>(const std::vector<std::string> &paths,
                                                 unsigned k, uint16_t w, spvec_t spaces,
                                                 void *data, int num_threads, size_t np);
template size_t estimate_cardinality<score::Hash>(const std::vector<std::string> &paths,
                                                  unsigned k, uint16_t w, spvec_t spaces,
                                                  void *data, int num_threads, size_t np);
template size_t estimate_cardinality<score::Entropy>(const std::vector<std::string> &paths,
                                                     unsigned k, uint16_t w, spvec_t spaces,
                                                     void *data, int num_threads, size_t np);
}
