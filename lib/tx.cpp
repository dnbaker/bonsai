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
    const uint64_t num_occ(name_map_->n_occupied);
    fwrite(&num_occ, sizeof(num_occ), 1, fp);
    fwrite(&n_syn_,  sizeof(n_syn_), 1, fp);
    fwrite(&ceil_,   sizeof(ceil_), 1, fp);
    size_t nwritten(0);
    fprintf(fp, "%s\n", name_.data());
    for(khiter_t ki(0); ki != kh_end(name_map_); ++ki)
        if(kh_exist(name_map_, ki))
            fprintf(fp, "%s\t%u\n", kh_key(name_map_, ki), kh_val(name_map_, ki)), ++nwritten;
    LOG_DEBUG("nw: %zu. no: %zu. nb: %zu.\n", nwritten, num_occ, name_map_->n_buckets);
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
        fgets(ts, sizeof(ts) - 1, fp);
        ki = kh_put(name, name_map_, ts, &khr);
        p = strchr(ts, '\t');
        kh_val(name_map_, ki) = atoi(p + 1);
        *p = '\0';
        kh_key(name_map_, ki) = strdup(ts);
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
    if(name_ != other.name_) return false;
    khiter_t ki, ki2;
    for(ki = 0; ki != kh_end(tax_map_); ++ki) {
        if(kh_exist(tax_map_, ki)) {
            if((ki2 = kh_get(p, other.tax_map_, kh_key(tax_map_, ki))) == kh_end(other.tax_map_)) return false;
            if(kh_val(tax_map_, ki) != kh_val(other.tax_map_, ki2)) return false;
        }
    }
    for(ki = 0; ki != kh_end(name_map_); ++ki)
        if(kh_exist(name_map_, ki))
            if((ki2 = kh_get(name, other.name_map_, kh_key(name_map_, ki))) == kh_end(other.name_map_) ||
                    kh_val(name_map_, ki) != kh_val(other.name_map_, ki2))
                return false;
    return true;
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

uint64_t vec_popcnt(std::string &vec) {
    uint64_t *arr((uint64_t *)&vec[0]), ret(0);
    for(size_t i(0), e(vec.size() >> 3); i < e; ++i)  ret += popcount(arr[i]);
    for(auto i(vec.data() + (vec.size() & (~0x7ul))), end(vec.size() + vec.data());
        i < end; ++i) ret += popcount(*i);
    return ret;
}

std::string rand_string(size_t n) {
    std::string ret;
    ret.reserve(n);
    static const char set[] = "abcdefghijklmnopqrstuvwxyz123456";
    while(ret.size() < n) ret += set[rand() & 31];
    assert(ret.size() == n);
    return ret;
}

} // namespace emp
