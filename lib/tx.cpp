#include "lib/tx.h"

namespace emp {

void Taxonomy::add_node_impl(const char *node_name, const unsigned node_id, const unsigned parent) {
    khint_t ki;
    int khr;
    if(kh_get(p, tax_map_, node_id) != kh_end(tax_map_)) LOG_EXIT("Indistinct node_id %u given. Abort!\n", node_id);
    ki = kh_put(name, name_map_, node_name, &khr);
    if(khr < 0) LOG_EXIT("Insufficient memory for adding to name hash table.\n");
    kh_val(name_map_, ki) = node_id;
}

void Taxonomy::add_node(const char *node_name, const unsigned parent) {
    unsigned new_id(1);
    assert(kh_get(p, tax_map_, 1) != kh_end(tax_map_)); // Need to have the root in the tree.
    while(kh_get(p, tax_map_, new_id) != kh_end(tax_map_))
        new_id = rand();
    add_node_impl(node_name, new_id, parent);
}

void Taxonomy::write(const char *fn) {
    khash_write<khash_t(p)>(tax_map_, fn);
    FILE *fp(fopen(fn, "a+b"));
    //fwrite(tax_map_->n_buckets
    for(khiter_t ki(0); ki != kh_end(name_map_); ++ki) {
        if(!kh_exist(name_map_, ki)) continue;
        for(const char *p(kh_key(name_map_, ki)); *p; fputc(*p++, fp));
        fputc('\0', fp);
        fwrite(&kh_val(name_map_, ki), 1, sizeof(kh_val(name_map_, ki)), fp);
    }
    fclose(fp);
}

} // namespace emp
