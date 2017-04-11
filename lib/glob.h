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

    tree_adjudicator_t(std::uint64_t orig, khash_t(p) *taxmap, int recalculate=1):
        full_taxmap_(taxmap), original_tax_count_(orig), max_el_(roundup64(orig)), recalculate_scores_(recalculate) {}

    void process_subtree(std::size_t i);

    template<typename... U>
    void emplace_subtree(U&&... Args) {
        subtrees_.emplace_back(std::forward<U>(Args)...);
        process_subtree(subtrees_.size() - 1);
        // Add potential_node_t's from each subtree using default scoring
    }
    int write_header(FILE *fp) {
        return std::fputs("#New node id\tParent\tChildren\n", fp);
    }
    int write_new_node(FILE *fp, const potential_node_t &node) {
        kstring_t ks{0, 0, 0};
        const std::vector<std::uint64_t> &v(*node.bits_);
        tax_t new_id(std::rand());
        khiter_t ki;
        int khr;
        while((ki = kh_get(p, full_taxmap_, new_id)) != kh_end(full_taxmap_)) new_id = std::rand();
        ki = kh_put(p, full_taxmap_, new_id, &khr);
        kh_val(full_taxmap_, ki) = subtrees_[node.subtree_index_].parent_;
        ksprintf(&ks, "%u\t%u\t", new_id, subtrees_[node.subtree_index_].parent_);
        for(std::size_t i(0); i < original_tax_count_; ++i) {
            if(v[i >> 6] & (1ull << i)) ksprintf(&ks, "%u,", subtrees_[node.subtree_index_].taxes_[i]);
        }
        ks.s[ks.l - 1] = '\n';
        const int ret(std::fputs(ks.s, fp));
        free(ks.s);
        return ret;
    }
    int simple_write(FILE *fp) {
        std::size_t n(0);
        int ret(0);
        for(const auto &node: nodes_) {
            if(n >= max_el_) break;
            ret += write_new_node(fp, node);
        }
        return ret;
    }
};

bool used_as_lca(const tax_t tax, const std::string &fld);

} // namespace emp



#endif
