#ifndef _GLOB_H_
#define _GLOB_H_
#include "lib/util.h"
#include "lib/tree_climber.h"
#include "lib/bitmap.h"

namespace emp {

struct potential_node_t {
    const unsigned                    subtree_index_;
    const std::vector<std::uint64_t> *bits_; // View, does not own.
    std::uint64_t                     score_;
    bool operator<(const potential_node_t &other) const {
        return std::tie(score_, bits_, subtree_index_) <
               std::tie(other.score_, other.bits_, other.subtree_index_);
    }
    bool operator>(const potential_node_t &other) const {
        return std::tie(score_, bits_, subtree_index_) >
               std::tie(other.score_, other.bits_, other.subtree_index_);
    }
    bool operator==(const potential_node_t &other) const {
        return std::tie(score_, bits_, subtree_index_) ==
               std::tie(other.score_, other.bits_, other.subtree_index_);
    }
    potential_node_t(const std::vector<std::uint64_t> *bits, std::uint64_t score, unsigned subtree_index):
            subtree_index_(subtree_index), bits_(bits), score_(score) {
        if(bits == nullptr) throw std::runtime_error("bitmap pointer is null.");
    }
};


class tree_glob_t {
    using tax_path_map_t = std::unordered_map<std::uint32_t, std::forward_list<std::string>>;
public:
    const tax_t                                parent_;
    count::Counter<std::vector<std::uint64_t>> counts_;
    adjmap_t                                   fwd_;
    adjmap_t                                   rev_;
    std::vector<tax_t>                         taxes_;
private:
    const std::string                          parent_path_;
    khash_t(all)                              *acceptable_;


    static constexpr const char *KMER_SUFFIX = ".kmers.bin";

public:
    std::string make_parent(const std::string &fld, const tax_t parent);
    tree_glob_t(khash_t(p) *tax, tax_t parent, const std::string &fld, const Spacer &sp, tax_path_map_t &tpm,
                const std::unordered_map<tax_t, std::vector<tax_t>> &invert, int num_threads=16);
    void add(const tax_t tax) {
        taxes_.push_back(tax);
#if !NDEBUG
        if(taxes_.size() % 100 == 0) LOG_DEBUG("ZOMG size of taxes is %zu\n", taxes_.size());
#endif
    }

    template<typename Container>
    void add_children(Container &&c) { // &&
        for(const auto tax: c) add(tax); // Add all the taxes.
    }
    template<typename Container>
    void add_children(Container &c)  { // &
        for(const auto tax: c) add(tax); // Add all the taxes.
    }
    template<typename It>
    void add_children(It first, It end) {
        while(first != end)    add(*first);
    }
    void get_taxes(const std::unordered_map<tax_t, std::vector<tax_t>> &invert);
    ~tree_glob_t();
};

struct tree_adjudicator_t {
    using node_gt = std::greater<potential_node_t>;

    std::vector<tree_glob_t>            subtrees_;
    std::set<potential_node_t, node_gt> nodes_;
    khash_t(p)                         *full_taxmap_;
    const std::uint64_t                 original_tax_count_; // Needed for scoring
    const std::uint64_t                 max_el_;             // Needed to know which nodes to add.
    int                                 recalculate_scores_; // If yes, recalculate scores
    kstring_t                           ks_;

    tree_adjudicator_t(std::uint64_t orig, khash_t(p) *taxmap, int recalculate=1, int doublings=1):
        full_taxmap_(taxmap), original_tax_count_(orig), max_el_(roundup64(orig) << doublings), recalculate_scores_(recalculate), ks_{0, 0, 0} {}

    void process_subtree(std::size_t i);

    template<typename... U>
    void emplace_subtree(U&&... Args) {
        subtrees_.emplace_back(std::forward<U>(Args)...);
        process_subtree(subtrees_.size() - 1);
        // Add potential_node_t's from each subtree using default scoring
    }
    int write_header(std::FILE *fp) {
        return std::fputs("#New node id\tParent\tChildren\n", fp);
    }
    int write_new_node(std::FILE *fp, const potential_node_t &node);
    void adjust_children(const std::vector<std::uint64_t> *bits, const unsigned subtree_index) {
        std::set<const std::vector<std::uint64_t> *> adjusted;
        adjust_children(bits, adjusted, subtree_index);
    }
    void adjust_children(const std::vector<std::uint64_t> *bits,
                         std::set<const std::vector<std::uint64_t> *> &adjusted,
                         const unsigned subtree_index) {
        if(adjusted.find(bits) != adjusted.end()) return;
        auto it(subtrees_[subtree_index].fwd_.find(bits));
        if(it != subtrees_[subtree_index].fwd_.end()) {
            for(auto ptr: it->second) {
                adjust_children(ptr, adjusted, subtree_index);
            }
        } else {
            // Somehow adjust the children's scores. The problem is that our adjacency map is from pointers to pointers, and it needs to be from nodes to nodes.
        }
    }
    int adjusting_write(std::FILE *fp);
    int simple_write(std::FILE *fp);

    ~tree_adjudicator_t() {free(ks_.s);}
};

bool used_as_lca(const tax_t tax, const std::string &fld);

} // namespace emp



#endif
