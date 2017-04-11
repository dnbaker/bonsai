#include "lib/glob.h"

namespace emp {

std::string tree_glob_t::make_parent(const std::string &fld, const tax_t parent) {
    return std::string(fld.empty() ? "": fld + '/') + std::to_string(parent) + KMER_SUFFIX;
}


void tree_glob_t::get_taxes(const std::unordered_map<tax_t, std::vector<tax_t>> &invert) {
    add_children(get_all_descendents(invert, parent_));
}


tree_glob_t::~tree_glob_t() {if(acceptable_) khash_destroy(acceptable_);}


tree_glob_t::tree_glob_t(khash_t(p) *tax, tax_t parent, const std::string &fld, const Spacer &sp, tax_path_map_t &tpm,
            const std::unordered_map<tax_t, std::vector<tax_t>> &invert, int num_threads):
    parent_(parent), parent_path_(make_parent(fld, parent_)),
    //acceptable_(tree::load_binary_kmerset(parent_path_.data()))
    acceptable_(nullptr)
{
    LOG_DEBUG("About to get tax ids for children.\n");
    if(acceptable_) LOG_DEBUG("About to make tax counter. Size of acceptable hash: %zu\n", kh_size(acceptable_));
    get_taxes(invert);
    LOG_DEBUG("Got taxes: %zu\n", taxes_.size());
    tax_path_map_t tmp;
    for(auto tax: taxes_) {
        const auto m(tpm.find(tax));
        if(m != tmp.end()) tmp.emplace(tax, m->second);
    }
    // Reorder taxes to correspond with bits in the bitmap.
    taxes_.clear();
    for(auto &pair: tmp) taxes_.push_back(pair.first);

    counts_ = std::move(bitmap_t(kgset_t(tmp, sp, num_threads, acceptable_)).to_counter());
    LOG_DEBUG("Size of counts: %zu\n", counts_.size());
    if(acceptable_) {
        khash_destroy(acceptable_), acceptable_ = nullptr;
    }
    fwd_ = std::move(adjmap_t(counts_, true));
    rev_ = std::move(adjmap_t(counts_, false));
}

void tree_adjudicator_t::process_subtree(std::size_t i) {
    indices_.emplace(subtrees_[i].parent_, i);
    tree_glob_t &g(subtrees_[i]);
    const std::size_t tax_size(g.taxes_.size());
    for(auto &pair: g.counts_) {
        const std::uint64_t bitdiff(tax_size - popcnt::vec_popcnt(pair.first));
        std::uint64_t ret(pair.second * bitdiff);
        const auto m(g.fwd_.find(pair.first));
        if(m == g.fwd_.end()) {
            nodes_.emplace(g.parent_, &pair.first, ret);
            LOG_WARNING("Node not found in adjmap. Maybe leaf? Continuing....\n");
            continue;
        }
        for(const auto i: m->second) ret += bitdiff * g.counts_.find(*i)->second;
        nodes_.emplace(g.parent_, &pair.first, ret);
    }
}

bool used_as_lca(const tax_t tax, const std::string &fld) {
    std::string path(std::string(fld.empty() ? "": fld + '/') + std::to_string(tax) + ".kmers.bin");
    LOG_DEBUG("Checking for lca binary file at %s\n", path.data());
    return isfile(path.data());
}

} //namespace emp
