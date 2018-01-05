#include "kgset.h"
namespace emp {
void kg_helper(void *data_, long index, int tid) {
    kg_data *data((kg_data *)data_);
    Encoder<lex_score> enc(data->sp_);
    khash_t(all) *h(data->core_[index]);
    gzFile fp(gzopen(data->paths_[index].data(), "rb"));
    if(!fp) LOG_EXIT("Could not open file at %s\n", data->paths_[index].data());
    kseq_t *ks(kseq_init(fp));
    u64 min;
    int khr;
    while(kseq_read(ks) >= 0) {
        enc.assign(ks);
        while(enc.has_next_kmer())
            if((min = enc.next_minimizer()) != BF)
                if(!data->acceptable_ || (kh_get(all, data->acceptable_, min) != kh_end(data->acceptable_)))
                    kh_put(all, h, min, &khr);
    }
    kseq_destroy(ks);
    gzclose(fp);
}

void kg_list_helper(void *data_, long index, int tid) {
    kg_list_data &data(*(kg_list_data *)data_);
    auto &list(*data.fl_[index]);
    khash_t(all) *h(data.core_[index]);
    Encoder<lex_score> enc(data.sp_);
    LOG_DEBUG("Size of list: %zu. Performing for index %ld of %zu\n", size(list), index, data.core_.size());
    for(const auto &path: list) {
        gzFile fp(gzopen(path.data(), "rb"));
        if(!fp) LOG_EXIT("Could not open file at %s\n", path.data());
        kseq_t *ks(kseq_init(fp));
        u64 min;
        int khr;
        while(kseq_read(ks) >= 0) {
            enc.assign(ks);
            while(enc.has_next_kmer()) {
                if((min = enc.next_minimizer()) != BF) {
#if !NDEBUG
                    if(!data.acceptable_ ||
                       (kh_get(all, data.acceptable_, min) != kh_end(data.acceptable_) &&
                        LOG_DEBUG("Found %s in acceptable hash. New size of hash: %zu\n", data.sp_.to_string(min).data(), kh_size(h))))
#else
                    if(!data.acceptable_ ||
                       (kh_get(all, data.acceptable_, min) != kh_end(data.acceptable_)))
#endif
                        kh_put(all, h, min, &khr);
                }
            }
        }
        kseq_destroy(ks);
        gzclose(fp);
    }
    
}

} // namespace emp
