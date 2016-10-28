#include "db.h"

namespace kpg {

template<>
void build_minimized_database<tax_score>(khash_t(64) *td_map, const Spacer &sp, std::vector<std::string> &paths, chm_t &ret) {
    ret.clear();
    Encoder<tax_score> enc(0, 0, sp, td_map);
    uint64_t kmer;
    for(auto &path: paths) {
        gzFile fp(gzopen(path.data(), "rb"));
        kseq_t *ks(kseq_init(fp));
        khint_t ki;
        if(kseq_read(ks) >= 0) {
            enc.assign(ks);
            while(enc.has_next_kmer()) {
                if((kmer = enc.next_minimizer()) != BF) {
                    if(unlikely((ki = kh_get(64, td_map, kmer)) == kh_end(td_map)))
                        LOG_EXIT("kmer not in taxdepth hash. This should never happen.\n");
                    ret[kmer] = (uint32_t)kh_val(td_map, ki);
                }
            }
                    
        }
        gzclose(fp);
        kseq_destroy(ks);
    }
}

template<>
void build_minimized_database<hash_score>(khash_t(64) *td_map, const Spacer &sp, std::vector<std::string> &paths, chm_t &ret) {
    ret.clear();
    Encoder<hash_score> enc(0, 0, sp, td_map);
    uint64_t kmer;
    for(auto &path: paths) {
        gzFile fp(gzopen(path.data(), "rb"));
        kseq_t *ks(kseq_init(fp));
        khint_t ki;
        if(kseq_read(ks) >= 0) {
            enc.assign(ks);
            while(enc.has_next_kmer()) {
                if((kmer = enc.next_minimizer()) != BF) {
                    if(unlikely((ki = kh_get(64, td_map, kmer)) == kh_end(td_map)))
                        LOG_EXIT("kmer not in taxdepth hash. This should never happen.\n");
                    ret[kmer] = (uint32_t)kh_val(td_map, ki);
                }
            }
                    
        }
        gzclose(fp);
        kseq_destroy(ks);
    }
}
    
std::map<int, void (*)(khash_t(64) *, const Spacer &, std::vector<std::string> &, chm_t &)> dbbuild_map{
   {score_scheme::LEX, build_minimized_database<lex_score>},
   {score_scheme::TAX_DEPTH, build_minimized_database<hash_score>}
};

}
