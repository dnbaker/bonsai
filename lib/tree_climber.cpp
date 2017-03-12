#include "lib/tree_climber.h"

namespace emp { namespace tree {

void mw_helper(void *data_, long index, int tid);

struct safe_fp_t;
struct multi_writer_t {
    struct safe_fp_t {
        FILE *fp_;
        std::uint64_t n_;
        std::unique_ptr<std::mutex> m_;
        std::uint64_t a_[4096];
        int l_;
        safe_fp_t(const char *fname): fp_(fopen(fname, "wb")), n_(0), m_(std::move(new std::mutex())) {
            if(!fp_) LOG_EXIT("Could not open fname %s\n", fname);
            LOG_DEBUG("Opened file at %s\n", fname);
        }
        int write() {
            LOG_DEBUG("Writing %i elements\n", l_);
            const int ret(l_ ? std::fwrite(a_, l_, 8, fp_): 0);
            l_ = 0;
            return ret;
        }
        ~safe_fp_t() {
            write();
            fclose(fp_);
        }
        void push(std::uint64_t val) {
            LOG_DEBUG("Pushing back %" PRIu64 "\n", val);
            std::unique_lock<std::mutex> lock(*m_);
            a_[l_++] = val;
            if(l_ >= 4096) write();
        }
        safe_fp_t(safe_fp_t &&other): fp_(other.fp_), n_(other.n_), m_(std::move(other.m_)) {}
    };
    std::shared_mutex                    m_; // Hash insertion lock.
    std::unordered_map<tax_t, safe_fp_t> map_;
    std::vector<std::string>             paths_;
    const std::string                   &fld_;
    khash_t(c)                          *h_;
    int                                  num_threads_;
    std::size_t                          chunk_size_;
    std::atomic<std::uint64_t>           n_;
    auto find_or_insert(tax_t tax) {
        std::shared_lock<std::shared_mutex> lock(m_);
        auto m(map_.find(tax));
        if(m == map_.end()) {
            LOG_DEBUG("New element\n");
            std::unique_lock<std::shared_mutex> lock(m_);
            char buf[256];
            std::sprintf(buf, "%s%u.kmers.bin", fld_.data(), tax);
            paths_.emplace_back(buf);
            m = map_.emplace(tax, buf).first;
        } else LOG_DEBUG("STUFF\n");
        return m;
    }
    void add_entry(std::uint64_t index) {
        LOG_DEBUG("Adding entry\n");
        if(!kh_exist(h_, index)) return;
        auto iter(find_or_insert(kh_key(h_, index)));
        iter->second.push(h_->keys[index]);
        ++n_;
        LOG_DEBUG("Added entry\n");
    }
    multi_writer_t(khash_t(c) *h, const std::string &fld, int num_threads=8, std::size_t chunk_size=1<<16):
        fld_(fld), h_(h), num_threads_(num_threads), chunk_size_(chunk_size), n_(0) {}
};
    struct mw_helper_t {
        multi_writer_t &mw_;
        mw_helper_t(multi_writer_t &mw, std::size_t chunk_size=1<<16): mw_(mw) {}
    };

void mw_helper(void *data_, long index, int tid) {
    mw_helper_t &mh(*(mw_helper_t *)data_);
    LOG_INFO("launch helper %ld of %u and total to do is %zu\n", index, tid, kh_size(mh.mw_.h_) / mh.mw_.chunk_size_);
    for(std::uint64_t i(index * mh.mw_.chunk_size_), e(std::min(kh_end(mh.mw_.h_), i + mh.mw_.chunk_size_));
        i < e; ++i) {
        LOG_DEBUG("About to add entry\n");
        mh.mw_.add_entry(index);
        LOG_DEBUG("Added entry at %zu in range %zu to %zu\n", i, e, e - mh.mw_.chunk_size_);
    }
}
void add_all_entries(multi_writer_t &mw) {
    LOG_INFO("make helper\n");
    mw_helper_t helper(mw, mw.chunk_size_);
    LOG_INFO("made helper\n");
    kt_for(mw.num_threads_, mw_helper, (void *)&helper, (kh_size(mw.h_) + mw.chunk_size_ - 1) / mw.chunk_size_);
    LOG_INFO("All done!\n");
}

std::vector<std::string> par_invert(Database<khash_t(c)> &db, const char *folder, int num_threads, std::size_t chunk_size) {
    multi_writer_t mw(db.db_, folder, num_threads, chunk_size);
    add_all_entries(mw);
    return std::move(mw.paths_);
}


std::vector<std::string> invert_lca_map(Database<khash_t(c)> &db, const char *folder, int prebuilt) {
    std::vector<std::string> ret;
    std::unordered_map<tax_t, std::FILE *> ofps;
    std::unordered_map<tax_t, std::uint64_t> ofp_counts;
    khash_t(c) *map(db.db_); // Does not own; shorthand.
    assert(map);
    std::string fld;
    int fwritten;
#if !NDEBUG
    std::size_t n_processed(0);
#endif
    if(folder && *folder) fld = folder, fld += '/';
    if(prebuilt == 0) {
        for(khiter_t ki(0); ki != kh_end(map); ++ki) {
            if(!kh_exist(map, ki)) continue;
            auto m(ofps.find(kh_val(map, ki)));
            auto mc(ofp_counts.find(kh_val(map, ki)));
            if(m == ofps.end()) {
                LOG_DEBUG("Adding %u to map\n", kh_val(map, ki));
                const std::uint64_t tmp(UINT64_C(-1));
                char buf[256];
                std::sprintf(buf, "%s%u.kmers.bin", fld.data(), kh_val(map, ki));
                ret.emplace_back(buf);
                m = ofps.emplace(kh_val(map, ki), std::fopen(buf, "wb")).first;
                if(m->second == nullptr) {
                    char buf2[512];
                    std::sprintf(buf2, "Could not open file at path %s\n", buf);
                    perror(buf2);
                    exit(1);
                }
                LOG_DEBUG("Opened file at %s\n", buf);
                mc = ofp_counts.emplace(kh_val(map, ki), 1).first;
                fwritten = std::fwrite(&tmp, sizeof(tmp), 1, m->second);
                if(fwritten != 1) {
                    char buf2[256];
                    std::sprintf(buf2, "Could not write dummy value. Size written: %i\n", fwritten);
                    perror(buf), exit(1);
                }
            } else ++mc->second;
            fwritten = std::fwrite(&kh_key(map, ki), sizeof(kh_key(map, ki)), 1, m->second);
            if(fwritten != 1) {
                char buf2[256];
                std::sprintf(buf2, "Could not write dummy value. Size written: %i. taxid: %u\n", fwritten, kh_val(map, ki));
                perror(buf2), exit(1);
            }
#if !NDEBUG
            if((++n_processed & 0xFFFFFF) == 0) LOG_DEBUG("Processed %zu so far\n", n_processed);
#endif
        }
        for(auto &pair: ofps) {
            std::rewind(pair.second);
            std::fwrite(&ofp_counts[pair.first], 1, sizeof(std::uint64_t), pair.second);
            std::fclose(pair.second);
        }
    } else {
        LOG_DEBUG("Just getting filenames\n");
        char buf[256];
        for(khiter_t ki(0); ki != kh_end(map); ++ki) {
            if(!kh_exist(map, ki)) continue;
            auto m(ofps.find(kh_val(map, ki)));
            if(m == ofps.end()) {
                std::sprintf(buf, "%s%u.kmers.bin", fld.data(), kh_val(map, ki));
                ret.emplace_back(buf);
            }
        }
    }
    LOG_DEBUG("Got 'em!\n");
    return ret;
}

khash_t(p) *pruned_taxmap(std::vector<std::string> &paths, khash_t(p) *taxmap, khash_t(name) *name_hash) {
    std::set<tax_t> found, keep;
    std::uint64_t ki1, ki2;
    for(auto &path: paths) found.insert(get_taxid(path.data(), name_hash));
    for(auto i: found) {
        keep.insert(i);
        ki1 = kh_get(p, taxmap, i);
        i = kh_val(taxmap, ki1);
        while(i) {
            keep.insert(i);
            ki1 = kh_get(p, taxmap, i);
            i = kh_val(taxmap, ki1);
        }
    }
    int khr;
    khash_t(p) *ret(kh_init(p));
    for(const auto i: keep) {
        ki1 = kh_put(p, ret, i, &khr);
        ki2 = kh_get(p, taxmap, i);
        kh_val(ret, ki1) = kh_val(taxmap, ki2);
    }
    return ret;
}

khash_t(all) *load_binary_kmerset(const char *path) {
    std::FILE *fp(fopen(path, "rb"));
    khash_t(all) *ret(kh_init(all));
    std::uint64_t n;
    int khr;
    std::fread(&n, 1, sizeof(n), fp);
    kh_resize(all, ret, n);
    while(std::fread(&n, 1, sizeof(std::uint64_t), fp) == sizeof(std::uint64_t))
        kh_put(all, ret, n, &khr);
    std::fclose(fp);
    return ret;
}


std::vector<std::uint64_t> load_binary_kmers(const char *path) {
    std::FILE *fp(fopen(path, "rb"));
    std::vector<std::uint64_t> ret;
    std::uint64_t n, ind(0);
    std::fread(&n, sizeof(n), 1, fp);
    ret.resize(n); // Initializes to 0 unnecessarily. Better than passing it to a temporary variable every time, I think.
    while(std::fread(ret.data() + ind++, sizeof(std::uint64_t), 1, fp) == sizeof(std::uint64_t));
    std::fclose(fp);
    return ret;
}

} /* namespace tree */ } // namespace emp
