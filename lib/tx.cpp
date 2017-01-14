#include "lib/tx.h"

namespace emp {

void kg_helper(void *data_, long index, int tid) {
    kg_data *data((kg_data *)data_);
    Encoder<lex_score> enc(data->sp_);
    khash_t(all) *h(data->core_[index]);
    gzFile fp(gzopen(data->paths_[index].data(), "rb"));
    kseq_t *ks(kseq_init(fp));
    uint64_t min;
    int khr;
    while(kseq_read(ks) >= 0) {
        enc.assign(ks);
        while(enc.has_next_kmer())
            if((min = enc.next_minimizer()) != BF)
                kh_put(all, h, min, &khr);
    }
    kseq_destroy(ks);
    gzclose(fp);
}

void Taxonomy::add_node_impl(const char *node_name, const unsigned node_id, const unsigned parent) {
    khint_t ki;
    int khr;
    if(kh_get(p, tax_map_, node_id) != kh_end(tax_map_)) LOG_EXIT("Indistinct node_id %u given. Abort!\n", node_id);
    ki = kh_put(name, name_map_, node_name, &khr);
    if(khr < 0) LOG_EXIT("Insufficient memory for adding to name hash table.\n");
    kh_val(name_map_, ki) = node_id;
}

void Taxonomy::add_node(const char *node_name, const unsigned parent) {
    if(!has_capacity()) throw std::logic_error("Attempted to add nodes past capacity. Abort!");
    unsigned new_id(1);
    assert(kh_get(p, tax_map_, 1) != kh_end(tax_map_)); // Need to have the root in the tree.
    while(kh_get(p, tax_map_, new_id) != kh_end(tax_map_))
        new_id = rand();
    add_node_impl(node_name, new_id, parent);
}

void Taxonomy::write(const char *fn) const {
    FILE *fp(fopen(fn, "wb"));
    const uint64_t nb(tax_map_->n_buckets);
    fwrite(&nb,     sizeof(nb), 1, fp);
    fwrite(&n_syn_, sizeof(n_syn_), 1, fp);
    fwrite(&ceil_,  sizeof(ceil_), 1, fp);
    for(khiter_t ki(0); ki != kh_end(name_map_); ++ki) {
        if(!kh_exist(name_map_, ki)) continue;
        for(const char *p(kh_key(name_map_, ki)); *p; fputc(*p++, fp));
        fputc('\n', fp);
        fwrite(&kh_val(name_map_, ki), sizeof(kh_val(name_map_, ki)), 1, fp);
    }
    khash_write_impl(tax_map_, fp);
    fclose(fp);
}

Taxonomy::Taxonomy(const char *path, unsigned ceil): name_map_(kh_init(name)) {
    int khr;
    khiter_t ki;
    FILE *fp(fopen(path, "rb"));
    char ts[1 << 10];
    uint64_t n;
    fread(&n,       sizeof(n),      1, fp);
    fread(&n_syn_,  sizeof(n_syn_), 1, fp);
    fread(&ceil_,   sizeof(ceil_),  1, fp);
    kh_resize(name, name_map_, n);

    uint32_t tmp; 
    for(uint64_t i(0); i < n; ++i) {
        fgets(ts, sizeof(ts) - 1, fp);
        fread(&tmp, sizeof(tmp), 1, fp);
        ki = kh_put(name, name_map_, ts, &khr);
        kh_key(name_map_, ki) = strdup(ts);
        kh_val(name_map_, ki) = tmp;
    }

    tax_map_ = khash_load_impl<khash_t(p)>(fp);
    fclose(fp);
    ceil_ = ceil ? ceil: tax_map_->n_buckets << 1;
}

template<typename T>
unsigned popcount(T val) {
    return __builtin_popcount(val);
}

template<>
unsigned popcount(unsigned long long val) {
    return __builtin_popcountll(val);
}

template<>
unsigned popcount(unsigned long val) {
    return __builtin_popcountl(val);
}

uint64_t vec_popcnt(std::vector<uint64_t> &vec) {
    uint64_t ret(vec[0]);
    for(size_t i(1), end(vec.size()); i < end; ++i) ret += popcount(vec[i]);
    return ret;
}

} // namespace emp
