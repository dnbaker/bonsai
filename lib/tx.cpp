#include "lib/tx.h"

namespace emp {

void kg_helper(void *data_, long index, int tid) {
    kg_data *data((kg_data *)data_);
    Encoder<lex_score> enc(data->sp_);
    khash_t(all) *h(data->core_[index]);
    gzFile fp(gzopen(data->paths_[index].data(), "rb"));
    if(!fp) LOG_EXIT("Could not open file at %s\n", data->paths_[index].data());
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
    const uint64_t num_occ(name_map_->n_occupied);
    fwrite(&num_occ, sizeof(num_occ), 1, fp);
    fwrite(&n_syn_,  sizeof(n_syn_), 1, fp);
    fwrite(&ceil_,   sizeof(ceil_), 1, fp);
    size_t nwritten(0);
    fprintf(fp, "%s\n", name_.data());
    for(khiter_t ki(0); ki != kh_end(name_map_); ++ki) {
        if(kh_exist(name_map_, ki)) {
            fprintf(fp, "%s\t%u\n", kh_key(name_map_, ki), kh_val(name_map_, ki)), ++nwritten;
            assert(kh_get(name, name_map_, kh_key(name_map_, ki)) != kh_end(name_map_));
            assert(kh_key(name_map_, ki));
        }
    }
    LOG_DEBUG("nw: %zu. no: %zu. cl: %zu. nb: %zu.\n", nwritten, num_occ, ceil_, name_map_->n_buckets);
    assert(nwritten == num_occ);
    khash_write_impl(tax_map_, fp);
    fclose(fp);
}

Taxonomy::Taxonomy(const char *path, unsigned ceil): name_map_(kh_init(name)) {
    int khr;
    khiter_t ki;
    FILE *fp(fopen(path, "rb"));
    char ts[1 << 10], *p;
    uint64_t n;
    fread(&n,       sizeof(n),      1, fp);
    fread(&n_syn_,  sizeof(n_syn_), 1, fp);
    fread(&ceil_,   sizeof(ceil_),  1, fp);
    kh_resize(name, name_map_, n);
    fgets(ts, sizeof(ts) - 1, fp);
    *(strchr(ts, '\n')) = '\0';
    name_ = ts;
    LOG_DEBUG("n: %zu. syn: %zu. ceil: %zu. name: %s\n", n, n_syn_, ceil_, name_.data());

    for(uint64_t i(0); i < n; ++i) {
        fgets(ts, sizeof(ts), fp);
        *(p = strchr(ts, '\t')) = '\0';
        ki = kh_put(name, name_map_, ts, &khr);
        kh_val(name_map_, ki) = atoi(p + 1);
        kh_key(name_map_, ki) = strdup(ts);
        assert(strcmp(kh_key(name_map_, ki), ts) == 0);
    }

    tax_map_ = khash_load_impl<khash_t(p)>(fp);
    fclose(fp);
    ceil_ = ceil ? ceil: tax_map_->n_buckets << 1;
}

bool Taxonomy::operator==(Taxonomy &other) const {
    if(n_syn_ != other.n_syn_) {LOG_DEBUG("syn\n"); return false;}
    if(ceil_ != other.ceil_) {LOG_DEBUG("ceil %zu, %zu\n", ceil_, other.ceil_); return false; }
    if(!_kh_eq(tax_map_, other.tax_map_)) return false;
    if(!_kh_eq(name_map_, other.name_map_)) return false;
    if(name_ != other.name_) {
        LOG_DEBUG("name %s != %s\n", name_.data(), other.name_.data());
        return false;
    }

    for(khiter_t ki = 0; ki != kh_end(other.name_map_); ++ki) {
        if(kh_exist(other.name_map_, ki)) {
            if(kh_get(name, name_map_, kh_key(other.name_map_, ki)) == kh_end(name_map_)) {
                LOG_DEBUG("key %s missing from other\n", kh_key(name_map_, ki));
                return false;
            }
            if(kh_val(other.name_map_, ki) != kh_val(name_map_, kh_get(name, name_map_, kh_key(other.name_map_, ki)))) {
                LOG_DEBUG("Value is different.\n");
                return false;
            }
        }
    }

    return true;
}

std::string rand_string(size_t n) {
    std::string ret;
    ret.reserve(n);
    static const char set[] = "abcdefghijklmnopqrstuvwxyz123456";
    while(ret.size() < n) ret += set[rand() & 31];
    assert(ret.size() == n);
    return ret;
}

uint64_t vec_popcnt(const char *p, const size_t l) {
    uint64_t *arr((uint64_t *)p), ret(0);
    const uint64_t nloops((l + 7ul) >> 3);
    for(size_t i(0); i < nloops; ++i)  ret += popcount(arr[i]);
    return ret;
}

uint64_t vec_popcnt(const std::string &vec) {
    return vec_popcnt(vec.data(), vec.size());
}

} // namespace emp
