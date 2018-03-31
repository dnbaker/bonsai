#include "kgset.h"
namespace bns {
void kg_helper(void *data_, long index, int tid) {
    kg_data *data((kg_data *)data_);
    khash_t(all) *hash(data->core_[index]);
    int khr;
    Encoder<score::Lex> enc(data->sp_, data->canon_);
    LOG_INFO("Getting kmers from %s with index %ld\n", data->paths_[index].data(), index);
    enc.for_each([&](u64 min) {
        //LOG_INFO("Kmer is %s\n", data->sp_.to_string(min).data());
        if(!data->acceptable_ || (kh_get(all, data->acceptable_, min) != kh_end(data->acceptable_)))
            kh_put(all, hash, min, &khr);
        assert(kh_size(hash));
#if !NDEBUG
        
#endif
    }, data->paths_[index].data());
    LOG_INFO("kg helper! for path %s and thread id %i, I now have %zu kmers loaded.\n", data->paths_[index].data(), tid, size_t(kh_size(hash)));
}

void kg_list_helper(void *data_, long index, int tid) {
    kg_list_data &data(*(kg_list_data *)data_);
    auto &list(*data.fl_[index]);
    LOG_INFO("Size of list: %zu. Performing for index %ld of %zu\n", size(list), index, data.core_.size());
    khash_t(all) *hash(data.core_[index]);
    int khr;
    Encoder<score::Lex> enc(data.sp_, data.canon_);
    enc.for_each([&](u64 min) {
        if(!data.acceptable_ || (kh_get(all, data.acceptable_, min) != kh_end(data.acceptable_)))
            khash_put(hash, min, &khr);
    }, list);
}

} // namespace bns
