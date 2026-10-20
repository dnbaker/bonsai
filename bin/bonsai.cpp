#include <fstream>
#include <numeric>
#include <sstream>
#include <omp.h>
#include "bonsai/feature_min.h"
#include "bonsai/util.h"
#include "bonsai/database.h"
#include "bonsai/classifier.h"
#include "bonsai/bitmap.h"
#include "bonsai/tx.h"
#include "bonsai/setcmp.h"
#include "bonsai/flextree.h"
#include "bonsai/grafting.h"
#include "bonsai/kraken_export.h"
#include "bonsai/profile_stats.h"
#include "sketch/hll.h"
// #include "sketch/filterhll.h"

#if USE_CPPITERTOOLS
#include "cppitertools/groupby.hpp"
#endif

using namespace bns;

namespace {
inline double now_s() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}
}

using namespace std::literals;
using std::cerr;
using std::cout;
using std::begin;
using std::end;
template<typename SketchType>
inline hll::hll_t &get_hll(SketchType &s);
template<>
inline hll::hll_t &get_hll<hll::hll_t>(hll::hll_t &s) {return s;}
/*
template<>
inline hll::hll_t &get_hll<fhll::fhll_t>(fhll::fhll_t &s) {return s.hll();}
template<>
inline hll::hll_t &get_hll<fhll::pcfhll_t>(fhll::pcfhll_t &s) {return s.hll();}
*/

template<typename SketchType> inline const hll::hll_t &get_hll(const SketchType &s);
template<>
inline const hll::hll_t &get_hll<hll::hll_t>(const hll::hll_t &s) {return s;}
/*
template<>
inline const hll::hll_t &get_hll<fhll::fhll_t>(const fhll::fhll_t &s) {return s.hll();}
*/

template<typename SketchType>
struct est_helper {
    const Spacer                      &sp_;
    const std::vector<std::string> &paths_;
    std::mutex                         &m_;
    const u64                          np_;
    const bool                      canon_;
    void                            *data_;
    std::vector<SketchType>         &hlls_;
    kseq_t                            *ks_;
};
template<typename SketchType, typename ScoreType=score::Lex>
void est_helper_fn(void *data_, long index, int tid) {
    est_helper<SketchType> &h(*(est_helper<SketchType> *)(data_));
    fill_lmers<ScoreType, SketchType>(h.hlls_[tid], h.paths_[index], h.sp_, h.canon_, h.data_, h.ks_ + tid);
}

template<typename SketchType, typename ScoreType=score::Lex>
void fill_sketch(SketchType &ret, const std::vector<std::string> &paths,
              unsigned k, uint16_t w, const spvec_t &spaces, bool canon=true,
              void *data=nullptr, int num_threads=1, u64 np=23, kseq_t *ks=nullptr) {
    // Default to using all available threads if num_threads is negative.
    if(num_threads < 0) {
        num_threads = std::thread::hardware_concurrency();
        LOG_INFO("Number of threads was negative and has been adjusted to all available threads (%i).\n", num_threads);
    }
    const Spacer space(k, w, spaces);
    if(num_threads <= 1) {
        for(u64 i(0); i < paths.size(); fill_lmers<ScoreType, SketchType>(ret, paths[i++], space, canon, data, ks));
    } else {
        std::mutex m;
        KSeqBufferHolder kseqs(num_threads);
        std::vector<SketchType> sketches;
        while(sketches.size() < (unsigned)num_threads) sketches.emplace_back(ret.clone());
        if(ret.size() != sketches.back().size()) throw "a party";
        est_helper<SketchType> helper{space, paths, m, np, canon, data, sketches, kseqs.data()};
        {
            ForPool pool(num_threads);
            pool.forpool(&est_helper_fn<SketchType, ScoreType>, &helper, paths.size());
        }
        auto &rhll = get_hll(ret);
        for(auto &sketch: sketches) rhll += get_hll(sketch);
    }
}

template<typename ScoreType=score::Lex>
hll::hll_t make_hll(const std::vector<std::string> &paths,
                unsigned k, uint16_t w, spvec_t spaces, bool canon=true,
                void *data=nullptr, int num_threads=1, u64 np=23, kseq_t *ks=nullptr, hll::EstimationMethod estim=hll::EstimationMethod::ERTL_MLE, uint16_t jestim=hll::JointEstimationMethod::ERTL_JOINT_MLE) {
    hll::hll_t global(np, estim, (hll::JointEstimationMethod)jestim);
    fill_sketch<hll::hll_t, ScoreType>(global, paths, k, w, spaces, canon, data, num_threads, np, ks);
    return global;
}

template<typename ScoreType=score::Lex>
u64 estimate_cardinality(const std::vector<std::string> &paths,
                            unsigned k, uint16_t w, spvec_t spaces, bool canon,
                            void *data=nullptr, int num_threads=-1, u64 np=23, kseq_t *ks=nullptr, hll::EstimationMethod estim=hll::EstimationMethod::ERTL_MLE) {
    auto tmp(make_hll<ScoreType>(paths, k, w, spaces, canon, data, num_threads, np, ks, estim));
    return tmp.report();
}

// bonsai compact <in.db> <out.cdb>: convert a khash classification database
// into the compact (cht_map_t) format; classify detects it automatically.
int compact_main(int argc, char *argv[]) {
    if(argc < 3) {
        std::fprintf(stderr, "Usage: %s <in.db> <out.cdb>\nConvert a classification database to the compact "
                             "(sketch cht_map_t) format. classify detects the format automatically.\n", *argv);
        return EXIT_FAILURE;
    }
    Database<khash_t(c)> db(argv[1]);
    // key domain (KmerKeys): 2k bits for 2-bit nucleotide k-mers (k <= 32),
    // 61 for rolling-hash keys (k > 32); checked against the keys present
    unsigned kb = KmerKeys<>::key_bits_for(db.k_);
    if(CompactKmerMap::key_bits_needed(db.db_) > kb) {
        LOG_WARNING("keys span more than the expected %u bits; using 64-bit keys\n", kb);
        kb = 64;
    }
    CompactKmerMap m(db.db_, CompactKmerMap::DEFAULT_MAX_LOAD, kb);
    m.finalize();
    CompactDatabase::write(argv[2], db.k_, db.w_, db.s_, m);
    LOG_INFO("Wrote compact database %s: %zu k-mers, %u-bit keys, %u-bit values, %zu bytes (khash: %zu bytes)\n",
             argv[2], m.size(), m.key_bits(), m.value_bits(), m.bytes_used(), khash_bytes(db.db_));
    return EXIT_SUCCESS;
}

bool endswith(const std::string &path, const std::string &suf) {
    return path.size() >= suf.size() && std::equal(std::crbegin(suf), std::crend(suf), std::crbegin(path));
}

int classify_main(int argc, char *argv[]) {
    int co, num_threads(1), emit_kraken(1), emit_fastq(0), emit_all(0), chunk_size(1 << 20), per_set(32);
    double confidence(0.), confidence_hi(0.);
    unsigned min_hit_groups(0);
    bool canonicalize(true), emit_evidence(false), hits_only(false), quality(false), placement(false);
    double species_identity(0.95), family_identity(0.85);
    bool uniq_p_set(false);
    std::string uniq_path, sketch_dir, profile_path;
    bool minhash(false);
    unsigned uniq_hll_p(12), sketch_b(16);
    std::ios_base::sync_with_stdio(false);
    std::FILE *ofp(stdout);
    if(argc < 4) {
        usage:
        std::fprintf(stderr, "Usage:\n%s <dbpath> <tax_path> <inr1.fq> [Optional: <inr2.fq>]\n"
                             "Flags:\n-o:\tRedirect output to path instead of stdout.\n"
                             "-c:\tSet chunk size. Default: %i\n"
                             "-a:\tEmit all records, not just classified.\n"
                             "-p:\tSet number of threads. [1] (Set -1 to use all threads.)\n"
                             "-k:\tEmit kraken-style output.\n"
                             "-K:\tDo not emit kraken-style output.\n"
                             "-f:\tEmit fastq-style output.\n"
                             "-x:\tConfidence threshold in [0, 1] (Kraken2's --confidence): walk the call up the\n"
                             "   \ttaxonomy until its clade holds this fraction of the read's non-ambiguous k-mers. [0]\n"
                             "-g:\tMinimum number of distinct hitting keys (hit groups) to classify a read. [0]\n"
                             "-X:\tSecond, high-confidence call (H: field): the -x walk-up of the call at this stricter\n"
                             "   \tthreshold in [0, 1]. [0: off]\n"
                             "-R:\tConfidence (-x, -X) over the read's hits only, not hits + misses.\n"
                             "-Q:\tConfidence (-x, -X) over the expected number of error-free k-mers given the read's\n"
                             "   \tbase qualities, not hits + misses: noisy reads from genomes in the database keep\n"
                             "   \ttheir call, distant relatives (hits only in conserved genes) do not. FASTQ only.\n"
                             "-I:\tIdentity-capped placement, \"species,family\" identity thresholds (e.g. 0.95,0.85): the read's\n"
                             "   \tidentity to the database is estimated from hits over expected error-free k-mers; below\n"
                             "   \tthe first the call is at most a family, below the second at most an order, so reads of\n"
                             "   \ttaxa missing from the database land at a supported rank. Needs ranks (nodes.dmp).\n"
                             "-U:\tAlso write a KrakenUniq-style report to this path: per taxon and clade, reads,\n"
                             "   \tk-mer hits, distinct k-mers (exact, then HyperLogLog), duplication and coverage of the\n"
                             "   \tdatabase's k-mers. Classification output is unchanged; off by default.\n"
                             "-u:\tSketch precision (2^p registers) for -U / -B. [12 HyperLogLog, 10 b-bit MinHash]\n"
                             "-B:\tWith -U: count with b-bit MinHash instead of HyperLogLog and write each clade's read\n"
                             "   \tk-mer sketch to this directory as <taxid>.bbmh (compare with reference sketches by\n"
                             "   \tcontainment / Jaccard).\n"
                             "-b:\tBits per register of the written -B sketches. [16]\n"
                             "-P:\tAlso write this sample's profile to this path (.bprof): per taxon its reads, k-mer hits\n"
                             "   \tand distinct k-mer counter, before clade folding, so profiles of read subsets merge\n"
                             "   \texactly and samples can be compared without the database (bonsai profile).\n"
                             "-m:\tCount distinct k-mers with b-bit MinHash (as -B does) rather than HyperLogLog: much\n"
                             "   \ttighter shared-k-mer estimates in bonsai profile compare.\n"
                             "-E:\tEmit lineage evidence of the call: E:exact,ancestor,descendant,conflict k-mer hits\n"
                             "   \tand I:estimated identity.\n"
                             "-K:\tDo not emit fastq-formatted output.\n"
                             "\nIf -f and -k are set, full kraken output will be contained in the fastq comment field."
                             "\n  Default: kraken-style only output.\n",
                 *argv, 1 << 14);
        std::exit(EXIT_FAILURE);
    }
    while((co = getopt(argc, argv, "Cc:p:o:S:x:g:X:I:U:u:B:b:P:ERQafFkKmh?")) >= 0) {
        switch(co) {
            case 'h': case '?': goto usage;
            case 'C': canonicalize = false; break;
            case 'a': emit_all = 1; break;
            case 'c': chunk_size = std::atoi(optarg); break;
            case 'F': emit_fastq  = 0; break;
            case 'f': emit_fastq  = 1; break;
            case 'K': emit_kraken = 0; break;
            case 'k': emit_kraken = 1; break;
            case 'p': num_threads = std::atoi(optarg); break;
            case 'o': ofp = std::fopen(optarg, "w"); break;
            case 'S': per_set = std::atoi(optarg); break;
            case 'x': confidence = std::atof(optarg); break;
            case 'g': min_hit_groups = unsigned(std::atoi(optarg)); break;
            case 'X': confidence_hi = std::atof(optarg); break;
            case 'E': emit_evidence = true; break;
            case 'R': hits_only = true; break;
            case 'Q': quality = true; break;
            case 'U': uniq_path = optarg; break;
            case 'u': uniq_hll_p = unsigned(std::atoi(optarg)); uniq_p_set = true; break;
            case 'B': sketch_dir = optarg; break;
            case 'b': sketch_b = unsigned(std::atoi(optarg)); break;
            case 'P': profile_path = optarg; break;
            case 'm': minhash = true; break;
            case 'I': placement = true;
                      if(std::sscanf(optarg, "%lf,%lf", &species_identity, &family_identity) != 2)
                          RUNTIME_ERROR("-I takes species,family identity thresholds, e.g. 0.95,0.85");
                      break;
        }
    }
    LOG_ASSERT(ofp);
    if(sketch_dir.size() && uniq_path.empty()) RUNTIME_ERROR("-B needs -U <report>");
    if(sketch_dir.size() && (sketch_b < 1 || sketch_b > 64)) RUNTIME_ERROR("-b must be in [1, 64]");
    switch(argc - optind) {
        default: goto usage;
        case 3:  LOG_DEBUG("Processing in single-end mode.\n"); break;
        case 4:  LOG_DEBUG("Processing in paired-end mode.\n"); break;
    }
    khash_t(p) *taxmap(build_parent_map(argv[optind + 1]));
    const rank_map_t ranks = placement ? build_rank_map(argv[optind + 1]): rank_map_t{};
    if(placement && ranks.empty()) RUNTIME_ERROR(std::string("-I needs ranks (nodes.dmp third field) in ") + argv[optind + 1]);
    const rank_map_t *rp = placement ? &ranks: nullptr;
    const bool collect = uniq_path.size() || profile_path.size();
    const rank_map_t report_ranks = collect ? (placement ? ranks: build_rank_map(argv[optind + 1])): rank_map_t{};
    auto rank_of = [&](tax_t t) -> uint8_t {
        const auto it = report_ranks.find(t);
        return it == report_ranks.end() ? uint8_t(RANK_UNKNOWN): it->second;
    };
    auto rank_name = [&](tax_t t) -> std::string {return report_ranks.empty() ? "-": profile_rank_name(rank_of(t));};
    // Parent lookups for the -U report and -P profile (sketch cht table).
    const CompactTaxMap ctaxmap = collect ? CompactTaxMap(taxmap): CompactTaxMap(1);
    unsigned db_k = 0;
    // Configure, classify, and (with -U) write the KrakenUniq-style report;
    // for_each_db_taxon(f) calls f(taxon) once per database k-mer (coverage).
    auto run = [&](auto &c, auto &&for_each_db_taxon) {
        c.set_confidence(confidence); c.set_min_hit_groups(min_hit_groups);
        c.set_confidence_hi(confidence_hi); c.set_emit_evidence(emit_evidence); c.set_confidence_hits_only(hits_only);
        c.set_confidence_quality(quality); c.set_identity_placement(rp, species_identity, family_identity);
        const bool bbmh = sketch_dir.size() || minhash;
        if(collect) c.enable_uniq(uniq_p_set ? uniq_hll_p: bbmh ? 10u: 12u, bbmh ? UniqSketch::BBMH: UniqSketch::HLL);
        // We can use optind + 3 for both single-end and paired-end mode since the argument at
        // index argc is null when argc - optind == 3.
        process_dataset(c, taxmap, argv[optind + 2], argv[optind + 3], ofp, chunk_size, per_set);
        std::unordered_map<tax_t, uint64_t> db_kmers;
        if(collect) for_each_db_taxon([&](tax_t t) {++db_kmers[t];});
        if(profile_path.size()) {   // before -U folds the counters
            Profile header;
            header.name = profile_path.substr(profile_path.find_last_of('/') + 1);
            if(endswith(header.name, ".bprof")) header.name.resize(header.name.size() - 6);
            header.k = db_k; header.canonical = canonicalize; header.has_ranks = !report_ranks.empty();
            header.fingerprint = db_fingerprint(db_kmers, db_k, canonicalize);
            write_profile(profile_path, header, *c.uniq(), ctaxmap, db_kmers, rank_of);
        }
        if(uniq_path.size()) {
            std::FILE *ufp = std::fopen(uniq_path.data(), "w");
            if(!ufp) RUNTIME_ERROR("could not open " + uniq_path);
            write_uniq_report(ufp, *c.uniq(), ctaxmap, &db_kmers, rank_name, sketch_dir.size() ? &sketch_dir: nullptr, sketch_b);
            std::fclose(ufp);
        }
    };
    if(CompactDatabase::is_compact(argv[optind])) {
        // Compact database (bonsai compact / build -Z): k-mer lookups batched per
        // read; a hashed-tier database's key width is its hash width.
        CompactDatabase db(argv[optind]);
        db_k = db.k_;
        if(db.wide()) {
            ClassifierGeneric<score::Lex, CompactKmerMap128> c(db.db128_.get(), db.s_, db.k_, db.k_, num_threads, emit_all, emit_fastq, emit_kraken, canonicalize, db.keyspec());
            run(c, [&](auto &&f) {db.db128_->for_each([&](auto, tax_t t) {f(t);});});
        } else {
            ClassifierGeneric<score::Lex, CompactKmerMap> c(db.db_.get(), db.s_, db.k_, db.k_, num_threads, emit_all, emit_fastq, emit_kraken, canonicalize, db.keyspec());
            run(c, [&](auto &&f) {db.db_->for_each([&](auto, tax_t t) {f(t);});});
        }
    } else {
        Database<khash_t(c)> db(argv[optind]);
        db_k = db.k_;
        ClassifierGeneric<score::Lex> c(db.db_, db.s_, db.k_, db.k_, num_threads,
                                       emit_all, emit_fastq, emit_kraken, canonicalize);
        run(c, [&](auto &&f) {
            for(khiter_t ki = 0; ki != kh_end(db.db_); ++ki) if(kh_exist(db.db_, ki)) f(kh_val(db.db_, ki));
        });
    }
    if(ofp != stdout) std::fclose(ofp);
    kh_destroy(p, taxmap);
    LOG_INFO("Successfully completed classify!\n");
    return EXIT_SUCCESS;
}


int phase2_main(int argc, char *argv[]) {
    int c, mode(score_scheme::LEX), wsz(-1), num_threads(1), k(31);
    bool canon(true), compact_out(false);
    KeySpec keyspec;              // -P alphabet, -6 translate, -H hash bits (compact builds)
    WRITE write_fmt = UNCOMPRESSED;
    std::size_t start_size(1<<16);
    std::string spacing, tax_path, seq2taxpath, paths_file;
    std::ios_base::sync_with_stdio(false);
    std::string dbpath;
    // TODO: update documentation for tax_path and seq2taxpath options.
    if(argc < 4) {
        usage:
        std::fprintf(stderr, "Usage: %s <flags> -T <tax.path> -M <seq2tax.path> <out.path> <paths>\n"
                     "(minimized modes instead take <phase1map.path> <out.path> <paths>)\nFlags:\n"
                     "-k: Set k.\n"
                     "-p: Number of threads [1] (set to -1 to use all threads)\n"
                     "-t: Build for taxonomic minimizing\n-f: Build for feature minimizing\n"
                     "-F: Load paths from file provided instead further arguments on the command-line.\n"
                     "-e: Use entropy maximization.\n"
                     "-f: Use feature count minimization.\n"
                     "-t: Use tax depth maximization.\n"
                     "-w: Set window size.\n"
                     "-T: Set tax_path.\n"
                     "-M: Set seq2taxpath.\n"
                     "-S: Set spacing.\n"
                     "-z: Write gzip-compressed.\n"
                     "-Z: Build a compact database (sketch cht_map_t; see `bonsai compact`) directly, without\n"
                     "    an intermediate khash. Lex LCA databases only; classify detects the format.\n"
                     "-P: Amino-acid key alphabet for a protein database (-Z): PROTEIN20, PROTEIN14, PROTEIN8, PROTEIN6.\n"
                     "    k counts residues. References are protein FASTA unless -6 is given.\n"
                     "-6: 6-frame translate nucleotide references (and, at classify, reads) into -P residues.\n"
                     "-H: Key width in bits for hashed-tier keys (-Z): 16..61 [61], or 65..128 for 128-bit keys.\n"
                     , *argv);
        std::exit(EXIT_FAILURE);
    }
    while((c = getopt(argc, argv, "Cw:M:S:p:k:T:F:P:H:tef6Zh?")) >= 0) {
        switch(c) {
            case 'C': canon = false; break;
            case 'h': case '?': goto usage;
            case 'k': k = std::atoi(optarg); break;
            case 'p': num_threads = std::atoi(optarg); break;
            case 'S': spacing = optarg; break;
            case 's': start_size = strtoull(optarg, nullptr, 10); break;
            case 't': mode = score_scheme::TAX_DEPTH; break;
            case 'f': mode = score_scheme::FEATURE_COUNT; break;
            case 'w': wsz = std::atoi(optarg); break;
            case 'T': tax_path = optarg; break;
            case 'M': seq2taxpath = optarg; break;
            case 'F': paths_file = optarg; break;
            case 'e': mode = score_scheme::ENTROPY; break;
            case 'z': write_fmt = ZLIB; break;
            case 'Z': compact_out = true; break;
            case 'P': if(!input_type_from_string(optarg, keyspec.rht)) LOG_EXIT("Unknown alphabet %s\n", optarg); break;
            case '6': keyspec.translate = true; break;
            case 'H': keyspec.hash_bits = unsigned(std::atoi(optarg)); break;
        }
    }
    dbpath = argv[optind];
    if(num_threads < 0) num_threads = std::thread::hardware_concurrency();
    if(compact_out && mode != score_scheme::LEX) LOG_EXIT("-Z (compact output) builds lex LCA databases only.\n");
    if((keyspec.rht != DNA || keyspec.translate || keyspec.hash_bits) && !compact_out)
        LOG_EXIT("-P / -6 / -H describe the key encoding, which only compact databases (-Z) record.\n");
    if(keyspec.translate && keyspec.rht == DNA) LOG_EXIT("-6 translates into an amino-acid alphabet; give one with -P.\n");
    if(keyspec.rht != DNA) canon = false;
    if(wsz < k) wsz = k;
    // k <= 32: exact 64-bit packed k-mers. 32 < k <= 64 (lex LCA databases):
    // 61-bit keys folded from the exact 128-bit k-mer (spaced seeds and windows
    // supported), or exact 128-bit keys with -Z. k > 64: rolling-hash keys,
    // unspaced only.
    // capacities are per alphabet (DNA: 32 / 64 per 64- / 128-bit word)
    if(k > int(rh2n(keyspec.rht, 8)) && mode != score_scheme::LEX)
        LOG_EXIT("k beyond a 64-bit key needs the default lex mode.\n");
    if(k > int(rh2n(keyspec.rht, 16)) && !spacing.empty())
        LOG_EXIT("k beyond a 128-bit key uses rolling-hash keys, which need an unspaced seed.\n");
#ifdef ZWRAP_USE_ZSTD
    const std::string suf(".zst");
#else
    const std::string suf(".gz");
#endif
    if(endswith(dbpath, suf))     write_fmt = ZLIB;
    if(write_fmt && !endswith(dbpath, ".gz"))
        dbpath += suf, LOG_INFO("Writing gzipped, but without a .gz suffix. Adding it.\n");
    LOG_INFO("db output path: %s\n", dbpath.data());
    spvec_t sv(parse_spacing(spacing.data(), k));
    // The lex/entropy branch takes positional <out.path> <paths...> (see the
    // README); the minimized branch below takes <phase1map> <out.path> <paths...>.
    std::vector<std::string> inpaths(paths_file.size() ? get_paths(paths_file.data())
                                                       : std::vector<std::string>(argv + optind + 1, argv + argc));
    if(inpaths.empty()) LOG_EXIT("Need input files from command line or file. See usage.\n");
    LOG_DEBUG("Got paths\n");
    if(seq2taxpath.empty()) LOG_EXIT("seq2taxpath required for final database generation.");
    if(score_scheme::LEX == mode || score_scheme::ENTROPY) {
        LOG_INFO("Final map will be written to %s\n", dbpath.data());
        Spacer sp(k, wsz, sv);
        Database<khash_t(c)>  phase2_map(sp);
        // Force using hll so that we can use __sync_bool_compare_and_swap to parallelize.
        LOG_INFO("About to estimate cardinality\n");
        // Cardinality of what will actually be stored: minimizers when the
        // window exceeds k, all k-mers otherwise. Estimating at w = k would
        // oversize the table by the minimization factor (e.g. 40x at
        // k = 25, w = 100 over real genomes).
        std::size_t hash_size(estimate_cardinality<score::Lex>(inpaths, k, wsz, sv, canon, nullptr, num_threads, 24));
#if !NDEBUG
        {
            uint64_t sum = 0;
            for(const auto &path: inpaths) {
                gzFile fp = gzopen(path.data(), "rb");
                kseq_t *ks = kseq_init(fp);
                while(kseq_read(ks) >= 0) {
                    sum += ks->seq.l;
                }
                kseq_destroy(ks);
                gzclose(fp);
            }
            assert(sum > hash_size || !std::fprintf(stderr, "sum: %" PRIu64". hash size: %zu\n", sum, hash_size));
        }
#endif
        LOG_INFO("Estimated cardinality: %zu\n", hash_size);
        if(tax_path.empty()) RUNTIME_ERROR("Tax path required. [See -T option.]");
        LOG_INFO("Parent map bulding from %s\n", tax_path.data());
        khash_t(p) *taxmap(build_parent_map(tax_path.data()));
        if(compact_out) {
            // key tiers: k <= 32 exact 64-bit; 32 < k <= 64 exact 128-bit (compact
            // databases only); k > 64 61-bit rolling-hash keys
            auto write_out = [&](const auto &m) {
                CompactDatabase::write(dbpath.data(), phase2_map.k_, phase2_map.w_, phase2_map.s_, *m, keyspec);
                LOG_INFO("Wrote compact database %s: %zu k-mers, %u-bit keys, %u-bit values, %zu bytes\n",
                         dbpath.data(), m->size(), m->key_bits(), m->value_bits(), m->bytes_used());
            };
            const unsigned rk = unsigned(k);
            if(rk > rh2n(keyspec.rht, 8) && rk <= rh2n(keyspec.rht, 16)) {   // exact 128-bit keys
                using K = KmerKeys<score::Lex, u128>;
                write_out(compact_lca_map<score::Lex, u128>(inpaths, taxmap, seq2taxpath.data(), sp, num_threads, canon, hash_size,
                          K::key_bits_for(rk, keyspec.rht, keyspec.hash_bits ? keyspec.hash_bits: K::DEFAULT_HASH_BITS), keyspec));
            } else {
                using K = KmerKeys<>;
                write_out(compact_lca_map<score::Lex>(inpaths, taxmap, seq2taxpath.data(), sp, num_threads, canon, hash_size,
                          K::key_bits_for(rk, keyspec.rht, keyspec.hash_bits ? keyspec.hash_bits: K::DEFAULT_HASH_BITS), keyspec));
            }
            kh_destroy(p, taxmap);
            return EXIT_SUCCESS;
        }
        //LOG_INFO("I just feel like stopping this executable now for testing.\n");
        //goto fail;
        phase2_map.db_ = score_scheme::LEX == mode ? lca_map<score::Lex>(inpaths, taxmap, seq2taxpath.data(), sp, num_threads, canon, hash_size)
                                                   : lca_map<score::Entropy>(inpaths, taxmap, seq2taxpath.data(), sp, num_threads, canon, hash_size);
        phase2_map.write(dbpath.data(), write_fmt);
        //fail:
        kh_destroy(p, taxmap);
        return EXIT_SUCCESS;
    }
    LOG_INFO("Making minimized map\n");
    // Minimized layout: <phase1map> <out.path> <paths...>; the first two
    // positionals are the input map and the output, not genomes.
    if(paths_file.empty()) inpaths.assign(argv + optind + 2, argv + argc);
    Database<khash_t(64)> phase1_map{Database<khash_t(64)>(dbpath.data())};
    Database<khash_t(c)>  phase2_map{phase1_map};
    Spacer sp(k, wsz, phase1_map.s_);
    khash_t(p) *taxmap(tax_path.empty() ? nullptr: build_parent_map(tax_path.data()));
    phase2_map.db_ = minimized_map<score::Hash>(inpaths, phase1_map.db_, seq2taxpath.data(), taxmap, sp, num_threads, start_size, canon);
    std::string dbpath2 = argv[optind + 1];
    if(endswith(dbpath2, suf))     write_fmt = ZLIB;
    if(write_fmt && !endswith(dbpath2, ".gz"))
        dbpath2 += suf, LOG_INFO("Writing gzipped, but without a .gz suffix. Adding it.\n");
    // Write minimized map
    phase2_map.write(dbpath2.data(), write_fmt);
    if(taxmap) kh_destroy(p, taxmap);
    return EXIT_SUCCESS;
}


int phase1_main(int argc, char *argv[]) {
    int c, taxmap_preparsed(0), use_hll(0), mode(score_scheme::LEX), wsz(-1), k(31), num_threads(1), sketch_size(24);
    bool canon(true);
    std::ios_base::sync_with_stdio(false);
    std::string spacing;

    if(argc < 5) {
        usage:
        std::fprintf(stderr, "Usage: %s <flags> <seq2tax.path> <taxmap.path> <out.path> <paths>\nFlags:\n"
                     "-k: Set k.\n"
                     "-p: Number of threads. [1] (Set to -1 to use all threads.)\n-S: add a spacer of the format "
                     "<int>,<int>,<int>, (...), where each integer is the number of spaces"
                     "between successive bases included in the seed. There must be precisely k - 1"
                     "elements in this list. Use this option multiple times to specify multiple seeds.\n"
                     "-s: add a spacer of the format <int>x<int>,<int>x<int>,"
                     "..., where the first integer corresponds to the space "
                     "between bases repeated the second integer number of times.\n"
                     "-S: Set HyperLogLog sketch size. For very large cardinalities, this may need to be increased for accuracy.\n"
                     "-t: Build for taxonomic minimizing.\n-f: Build for feature minimizing.\n"
                     "-H: Estimate rather than count kmers exactly before building map.\n"
                     "-T: Path to taxonomy map to load, if you've preparsed it. Not really worth it, building from scratch is fast.\n"
                     "-d: Write out in database format version 1.\n"
                     , *argv);
        std::exit(EXIT_FAILURE);
    }
    if("lca"s == argv[0])
        std::fprintf(stderr, "[W:%s] lca subcommand has been renamed phase1. "
                             "This has been deprecated and will be removed.\n", __func__);
    while((c = getopt(argc, argv, "Cs:S:p:k:tfTHh?")) >= 0) {
        switch(c) {
            case 'C': canon = false; break;
            case 'h': case '?': goto usage;
            case 'k': k = std::atoi(optarg); break;
            case 'p': num_threads = std::atoi(optarg); break;
            case 's': spacing = optarg; break;
            case 'S': sketch_size = std::atoi(optarg); break;
            case 'T': taxmap_preparsed = 1; break;
            case 'H': use_hll = 1; break;
            case 't': mode = score_scheme::TAX_DEPTH; break;
            case 'f': mode = score_scheme::FEATURE_COUNT; break;
            //case 'w': wsz = std::atoi(optarg); break;
        }
    }
    if(num_threads < 0) num_threads = std::thread::hardware_concurrency();
    wsz = k;
    if(k > 32) LOG_EXIT("k must be <= 32 for this alphabet (64-bit encodings).\n");
    khash_t(p) *taxmap(taxmap_preparsed ? khash_load<khash_t(p)>(argv[optind + 1])
                                        : build_parent_map(argv[optind + 1]));
    spvec_t sv(parse_spacing(spacing.data(), k));
    Spacer sp(k, wsz, sv);
    std::vector<std::string> inpaths(argv + optind + 3, argv + argc);
    std::size_t hash_size(use_hll ? estimate_cardinality<score::Lex>(inpaths, k, k, sv, canon, nullptr, num_threads, sketch_size): 1 << 16);
    if(use_hll) LOG_INFO("Estimated number of elements: %zu\n", hash_size);

    if(mode == score_scheme::LEX) LOG_EXIT("No phase1 required for lexicographic. Use phase2 instead.\n");
    auto mapbuilder(mode == score_scheme::TAX_DEPTH ? taxdepth_map<score::Lex>
                                                    : ftct_map<score::Lex>);
    Database<khash_t(64)> db(sp, 1, mapbuilder(inpaths, taxmap, argv[optind], sp, num_threads, canon, hash_size));
    for(auto &i: db.s_) {
        LOG_DEBUG("Decrementing value %i to %i\n", i, i - 1);
        --i;
    }
    db.write(argv[optind + 2]);

    kh_destroy(p, taxmap);
    return EXIT_SUCCESS;
}

int hist_main(int argc, char *argv[]) {
    for(const char *const *p(argv); *p; ++p) if((std::strcmp(*p, "-h") && std::strcmp(*p, "--help")) == 0) goto usage;
    if(argc == 1) {
        usage:
        std::fprintf(stderr, "Produces a histogram of how many kmers have been assigned to all taxids\n"
                             "Usage: bonsai %s <database.db> [outfile (omit to emit to stdout)]\n", *argv);
        std::exit(EXIT_FAILURE);
    }
    Database<khash_t(c)> db(argv[1]);
    khash_t(c) *map(db.db_);
    std::FILE *ofp(stdout);
    count::Counter<u32> counter;
    if(argc > 2) ofp = std::fopen(argv[2], "w");
    for(khiter_t ki(0); ki != kh_end(map); ++ki) if(kh_exist(map, ki)) counter.add(kh_val(map, ki));
    auto &cmap(counter.get_map());
    using elcount = std::pair<tax_t, u32>;
    std::vector<elcount> structs;
    for(auto& i: cmap) structs.emplace_back(i.first, i.second);
    SORT(std::begin(structs), std::end(structs), [] (const elcount &a, const elcount &b) {return a.second < b.second;});
    std::fputs("Name\tCount\n", ofp);
    for(const auto &i: structs) std::fprintf(ofp, "%u\t%u\n", i.first, i.second);
    if(ofp != stdout) std::fclose(ofp);
     return EXIT_SUCCESS;
 }

int err_main(int argc, char *argv[]) {
    std::fprintf(stderr, "[bonsai:%s] No valid subcommand provided. Options: prebuild/p1/phase, build/p2/phase2, classify, compact, metatree, profile\n", BONSAI_VERSION);
    return EXIT_FAILURE;
}

int zomg_main(int argc, char *argv[]) {
    // This is testing code which will eventually prepare the reformed taxonomy (without and with)
    // new nodes for downstream work and classificaiton.
    std::vector<std::string> paths {"1", "2", "3"};
    TaxonomyReformation tr("Foobar", paths, nullptr);
    return EXIT_FAILURE;
}

int metatree_usage(const char *arg) {
    std::fprintf(stderr, "Usage: %s <db.path> <taxmap> <nameidmap> <out_taxmap> <out_taxkey>\n"
                         "\n"
                         "-F: Parse file paths from file instead of further arguments at command-line.\n"
                         "-d: Do not perform inversion (assume it's already been done.)\n"
                         "-f: Store binary dumps in folder <arg>.\n"
                         "-L: Set an lca to restrict analysis to one portion of the subtree.\n"
                         "    Repeating this option multiple times accepts genomes which are descendents of any taxid provided by this option.\n"
                         "-p: nthreads [1] (set to -1 to use all threads.)\n"
                         "-n: nthreads [1]\n"
                         "-r: Reform the taxonomy so that every genome gets a node (TODO.md item 2), write\n"
                         "    <arg>.names.gz, <arg>.o2n.gz, and <arg>.n2o.gz, and emit placements against the reformed ids.\n"
                 , arg);
    std::exit(EXIT_FAILURE);
    return EXIT_FAILURE;
}

bool accepted_pass(const khash_t(p) *taxmap, const std::vector<tax_t> &accepted, tax_t id) {
    if(accepted.empty()) return true;
    for(const auto el: accepted) if(lca(taxmap, el, id) == el) return true;
    return false;
}

int metatree_main(int argc, char *argv[]) {
    if(argc < 5) metatree_usage(*argv);
    int c, num_threads(1), k(31), nelem(0);
    size_t heap_size = 1 << 15;
    std::vector<tax_t> accept_lcas;
    FILE *ofp(stdout);
    std::string paths_file, folder, spacing, reform_prefix;
    std::ios_base::sync_with_stdio(false);
    while((c = getopt(argc, argv, "L:p:w:k:s:f:F:r:n:o:h?")) >= 0) {
        switch(c) {
            case '?': case 'h':     return metatree_usage(*argv);
            case 'f': folder      = optarg;                       break;
            case 'k': k           = std::atoi(optarg);            break;
            case 'F': paths_file  = optarg;                       break;
            case 'o': ofp         = std::fopen(optarg, "w");      break;
            case 'p': num_threads = std::atoi(optarg);            break;
            case 'n': nelem       = std::strtoull(optarg, 0, 10); break;
            case 'L': accept_lcas.push_back(std::atoi(optarg));   break;
            case 'r': reform_prefix = optarg;                     break;
        }
    }
    if(num_threads <= 0) num_threads = std::thread::hardware_concurrency();
    if(k > 32) LOG_EXIT("k must be <= 32 for this alphabet (64-bit encodings).\n");
    Spacer sp(k, k, nullptr);
    omp_set_num_threads(num_threads);
    khash_t(name) *name_hash(build_name_hash(argv[optind + 2]));
    LOG_DEBUG("Parsed name hash.\n");
    LOG_DEBUG("Got inpaths. Now building parent map\n");
    khash_t(p) *taxmap(build_parent_map(argv[optind + 1]));
    LOG_DEBUG("Got taxmap. Now getting maxtax\n");
    tax_t max_tax(get_max_val(taxmap));
    LOG_DEBUG("Now getting tax depths\n");
    auto tax_depths(get_tax_depths(taxmap, argv[optind + 1]));
    LOG_DEBUG("Got tax depths\n");
    std::unordered_set<tax_t> used_taxes;
    std::vector<std::string> inpaths(paths_file.size() ? get_paths(paths_file.data())
                                                       : std::vector<std::string>(argv + optind + 5, argv + argc));
    std::unordered_set<std::string> save;
    for(size_t i(0); i < inpaths.size(); ++i) {
        const auto &path(inpaths[i]);
#if !NDEBUG
        if((i % 500) == 0) LOG_DEBUG("At index %zu/%zu, save size is %zu\n", i, inpaths.size(), save.size());
#endif
        tax_t id;
        if((id = get_taxid(path.data(), name_hash)) != UINT32_C(-1)) {
            if(accepted_pass(taxmap, accept_lcas, id)) {
                save.insert(path), used_taxes.insert(id);
            }
        }
    }
    {
        inpaths = std::vector<std::string>(save.begin(), save.end());
        std::unordered_set<std::string> tmp;
        std::swap(tmp, save);
    }
    std::vector<tax_t> raw_taxes(used_taxes.begin(), used_taxes.end());
    for(auto tax: raw_taxes) {
        tax_t id;
        while((id = get_parent(taxmap, tax)))
            used_taxes.insert(id), tax = id;
    }
    LOG_DEBUG("Finished filtering genomes\n");

    if(inpaths.size() == 0)
        throw std::runtime_error("No input paths. I need to process genomes to tell you about them.");
    LOG_INFO("Processing %zu genomes\n", inpaths.size());

// Core
    std::unordered_map<tax_t, std::forward_list<std::string>> tx2desc_map;
    std::vector<tax_t> taxes;
    std::unique_ptr<TaxonomyReformation> tr;
    if(!reform_prefix.empty()) {
        // Build a concise taxonomy (a node per genome) from the genomes,
        // their taxid annotations, and the NCBI taxonomy, then place k-mers
        // against the reformed ids.
        tr.reset(new TaxonomyReformation(argv[optind + 2], inpaths, taxmap));
        tr->write_name_map((reform_prefix + ".names.gz").data());
        tr->write_old_to_new((reform_prefix + ".o2n.gz").data());
        tr->write_new_to_old((reform_prefix + ".n2o.gz").data());
        taxmap = const_cast<khash_t(p) *>(tr->pmap());
        max_tax = get_max_val(taxmap);
        for(khiter_t ki(0); ki != kh_end(taxmap); ++ki)
            if(kh_exist(taxmap, ki) && kh_key(taxmap, ki) != 1) taxes.push_back(kh_key(taxmap, ki));
        for(const auto &pair: tr->newid_paths()) { // genome leaf -> path; walk ancestors
            for(tax_t t(pair.first), p(t); t; t = p) {
                auto it(tx2desc_map.find(t));
                if(it == tx2desc_map.end()) it = tx2desc_map.emplace(t, std::forward_list<std::string>{}).first;
                it->second.push_front(pair.second);
                const khint_t ki(kh_get(p, taxmap, t));
                p = ki == kh_end(taxmap) ? 0 : kh_val(taxmap, ki);
            }
        }
    } else {
        taxes = vector_set_filter(get_sorted_taxes(taxmap, argv[optind + 1]), used_taxes);
        tx2desc_map = tax2desc_genome_map(tax2genome_map(name_hash, inpaths), taxmap, taxes, tax_depths);
    }
#if !NDEBUG
    for(const auto tax: taxes) assert(kh_get(p, taxmap, tax) != kh_end(taxmap));
    ks::string ks;
    const int stderrfn(fileno(stderr));
    ks.resize(1 << 6);
    typename decltype(tx2desc_map)::iterator it;
    for(const auto tax: taxes) {
        ks.putsn("Tax ", 4);
        ks.putuw_(tax);
        ks.puts(" has for descendent genomes: ");
        if((it = tx2desc_map.find(tax)) == tx2desc_map.end()) ks.putsn_("N/A\n", 4);
        else {
            for(const auto el: it->second) {
                ks.putsn(el.data(), el.size()); ks.putc(',');
            }
            ks.back() = '\n';
        }
        if(ks.size() & (1 << 16)) ks.write(stderrfn), ks.clear();
    }
    ks.write(stderrfn);
    ks.clear();
#endif
    // (canonicalize, max_heapsz, to_add) -- positional drift here once bound
    // heap_size to canonicalize and nelem to max_heapsz_, crashing fill_heap.
    FMEmitter fme(taxmap, tx2desc_map, true, heap_size, nelem);
    // Group taxes by parent and emit per subtree; no cppitertools needed.
    pdqsort(taxes.begin(), taxes.end(), [taxmap](const tax_t a, const tax_t b) {
        return get_parent(taxmap, a) < get_parent(taxmap, b);
    });
    for(auto bit(taxes.begin()); bit != taxes.end();) {
        auto eit(bit);
        while(eit != taxes.end() && get_parent(taxmap, *eit) == get_parent(taxmap, *bit)) ++eit;
        fme.process_subtree(get_parent(taxmap, *bit), bit, eit, sp, num_threads, nullptr);
        bit = eit;
    }

    fme.run_collapse(max_tax, ofp);
    if(ofp != stdout) std::fclose(ofp);
    return EXIT_SUCCESS;
}


// First-token-of-header taxid lookup: reads the first fasta header and
// probes the name hash with its first whitespace-delimited token (a
// "|kraken:taxid|" suffix, if present, is stripped). Returns tax_t(-1) when
// the token is absent from the name hash.
// Create dir and any missing parents (mkdir -p semantics).
static void make_dirs(const std::string &dir) {
    std::string cur;
    cur.reserve(dir.size());
    for(const auto c: dir) {
        cur.push_back(c);
        if(c == '/' && cur.size() > 1) ::mkdir(cur.data(), 0755);
    }
    ::mkdir(dir.data(), 0755);
}

static tax_t export_header_taxid(const std::string &path, const khash_t(name) *name_hash, std::string &seqid_out) {
    const std::string seqid(first_header_token(path.data()));
    if(seqid.empty()) LOG_EXIT("Empty, unreadable, or headerless fasta file %s\n", path.data());
    seqid_out = seqid;
    const khint_t ki(kh_get(name, name_hash, seqid_out.data()));
    return ki == kh_end(name_hash) ? tax_t(-1): kh_val(name_hash, ki);
}

static std::unordered_map<tax_t, std::string> load_taxon_strings(const std::string &path) {
    std::unordered_map<tax_t, std::string> ret;
    std::FILE *fp(std::fopen(path.data(), "r"));
    if(fp == nullptr) LOG_EXIT("Could not open %s\n", path.data());
    char *line(nullptr);
    size_t cap(0);
    while(getline(&line, &cap, fp) >= 0) {
        char *tab(std::strchr(line, '\t'));
        if(tab == nullptr) continue;
        *tab = 0;
        char *end(tab + 1);
        while(*end && *end != '\n' && *end != '\r') ++end;
        *end = 0;
        ret[std::atoi(line)] = tab + 1;
    }
    std::free(line);
    std::fclose(fp);
    return ret;
}

int export_usage(const char *arg) {
    std::fprintf(stderr, "Usage: %s <flags> -T <taxmap> -M <names> -o <outdir> [genome.fa ...]\n"
                         "Emits nodes.dmp/names.dmp/seqid2taxid.map and kraken library directories\n"
                         "that kraken2-build and krakenuniq-build ingest directly (see\n"
                         "include/bonsai/kraken_export.h). With -B > 0, greedily grafts up to B\n"
                         "synthetic nodes (candidate lanes: observed k-mer classes + BBMinHash\n"
                         "tree-layered sketches) before exporting.\n"
                         "\n"
                         "-T: Taxonomy map (child->parent; NCBI nodes.dmp-style rows accepted).\n"
                         "-M: Name map: \"<seqid>\\t<taxid>\" per line (krakenuniq seqid2taxid.map format).\n"
                         "-F: Read genome paths from file instead of the command line.\n"
                         "-o: Output directory (required).\n"
                         "-k: k-mer size for the grafting fill. [%u]\n"
                         "-B: Graft budget (0 = export the tree as-is). [0]\n"
                         "-P: Synthetic node name prefix. [grafted_]\n"
                         "-N: Optional \"<taxid>\\t<name>\" file for taxonomy node names.\n"
                         "-R: Optional \"<taxid>\\t<rank>\" file for taxonomy node ranks.\n"
                         "\n"
                         "PlacementTree/KmerClasses use u64 genome masks: at most 64 genomes per export.\n",
                 arg, 25);
    std::exit(EXIT_FAILURE);
    return EXIT_FAILURE;
}

int export_main(int argc, char *argv[]) {
    int c, k(25), budget(0);
    std::string tax_path, name_path, outdir, paths_file, graft_prefix("grafted_"), names_file, ranks_file;
    std::ios_base::sync_with_stdio(false);
    while((c = getopt(argc, argv, "T:M:F:o:k:B:P:N:R:h?")) >= 0) {
        switch(c) {
            case 'h': case '?': return export_usage(*argv);
            case 'T': tax_path = optarg;       break;
            case 'M': name_path = optarg;      break;
            case 'F': paths_file = optarg;     break;
            case 'o': outdir = optarg;         break;
            case 'k': k = std::atoi(optarg);   break;
            case 'B': budget = std::atoi(optarg); break;
            case 'P': graft_prefix = optarg;   break;
            case 'N': names_file = optarg;     break;
            case 'R': ranks_file = optarg;     break;
        }
    }
    if(tax_path.empty() || name_path.empty() || outdir.empty()) return export_usage(*argv);
    if(k > 32) LOG_EXIT("k must be <= 32 for this alphabet (64-bit encodings).\n");
    std::vector<std::string> paths(paths_file.size() ? get_paths(paths_file.data())
                                                     : std::vector<std::string>(argv + optind, argv + argc));
    if(paths.empty()) return export_usage(*argv);
    khash_t(p) *taxmap(build_parent_map(tax_path.data()));
    khash_t(name) *name_hash(build_name_hash(name_path.data()));
    std::vector<std::string> seqids, kept;
    std::vector<tax_t> genome_taxa;
    for(const auto &path: paths) {
        std::string seqid;
        const tax_t tid(export_header_taxid(path, name_hash, seqid));
        if(tid == tax_t(-1)) {
            LOG_WARNING("Genome %s (seqid %s) is not in the name map; skipping.\n", path.data(), seqid.data());
            continue;
        }
        if(kh_get(p, taxmap, tid) == kh_end(taxmap)) {
            LOG_WARNING("Taxid %u for genome %s is not in the taxonomy; skipping.\n", tid, path.data());
            continue;
        }
        seqids.push_back(seqid);
        kept.push_back(path);
        genome_taxa.push_back(tid);
    }
    destroy_name_hash(name_hash);
    if(kept.empty()) LOG_EXIT("No genomes survived name/taxonomy filtering.\n");
    if(kept.size() > 64) LOG_EXIT("%zu genomes provided, but the placement tree is limited to 64 (u64 genome masks).\n", kept.size());
    paths.swap(kept);
    LOG_INFO("Exporting %zu genomes against %zu taxonomy nodes.\n", paths.size(), kh_size(taxmap));
    PlacementTree tree(taxmap, genome_taxa);
    if(budget > 0) {
        Spacer sp(k, k);
        KmerClasses kc;
        const double t0(now_s());
        kc.fill<score::Lex>(paths, sp);
        const double fill_t(now_s() - t0);
        const u64 delta0(tree.total_delta(kc));
        std::vector<u64> cands;
        for(const auto &pair: kc.count) if(pc64(pair.first) >= 2) cands.push_back(pair.first);
        const size_t nclasses(cands.size());
        const auto sks(genome_sketchers<score::Lex>(paths, sp));
        for(const auto &comp: sketch_cluster_candidates_tree(sks, tree, 0.008)) {
            u64 m(0);
            for(const auto g: comp) m |= u64(1) << g;
            cands.push_back(m);
        }
        LOG_INFO("Grafting: %zu classes, %zu sketch candidates, index fill %.3f s, initial delta %llu.\n",
                 nclasses, cands.size() - nclasses, fill_t, (unsigned long long)delta0);
        const GraftingResult gr(graft_greedy(tree, kc, cands, u32(budget)));
        const HammingDecomp hd(hamming_decomp(tree, kc));
        LOG_INFO("Grafted %zu nodes: delta %llu -> %llu (fn %llu, fp %llu).\n",
                 gr.added_masks.size(), (unsigned long long)delta0,
                 (unsigned long long)gr.deltas.back(), (unsigned long long)hd.fn, (unsigned long long)hd.fp);
    }
    std::unordered_map<tax_t, std::string> taxon_names, taxon_ranks;
    if(!names_file.empty()) taxon_names = load_taxon_strings(names_file);
    if(!ranks_file.empty()) taxon_ranks = load_taxon_strings(ranks_file);
    make_dirs(outdir);
    const KrakenExportResult result(KrakenExporter::write(tree, seqids, paths, outdir.data(), graft_prefix.data(),
                                                          taxon_names.empty() ? nullptr: &taxon_names,
                                                          taxon_ranks.empty() ? nullptr: &taxon_ranks));
    kh_destroy(p, taxmap);
    return result.seqid2taxid.empty() ? EXIT_FAILURE: EXIT_SUCCESS;
}

// bonsai profile: tools over classify -P sample profiles (see bonsai/profile.h).
namespace {
int profile_usage(const char *arg) {
    std::fprintf(stderr, "Usage: %s <subcommand> ...  (profiles: classify -P <out.bprof>)\n"
                 "  report [-B dir] [-b bits] <in.bprof> [out.tsv]\n"
                 "      The classify -U report (identical to classify's) recomputed from a profile; -B writes each\n"
                 "      clade's read k-mer sketch to dir/<taxid>.bbmh (b-bit MinHash profiles, classify -m or -B).\n"
                 "  merge [-n name] <out.bprof> <in.bprof>...\n"
                 "      Union of profiles of one database: lanes or read subsets of a sample (the result is the\n"
                 "      profile of all their reads, byte for byte), or a pooled group. Name: out's basename.\n"
                 "  compare [-r rank] [-m min_reads] [-K kmers.tsv] <in.bprof>...\n"
                 "      Every pair of samples, to stdout: Bray-Curtis and presence-Jaccard distances over taxa\n"
                 "      (clades at -r rank, e.g. species or genus; default: the taxa reads were called to) and the\n"
                 "      taxonomy earth mover's distance (unit edges). -K: per clade (at -r rank, default all) with\n"
                 "      >= min_reads reads in both samples, distinct read k-mers in each and shared (Jaccard,\n"
                 "      containments): same strain vs different strains of one species. -m: [1]\n"
                 "  table [-r rank] [-m min_reads] [-u] <in.bprof>...\n"
                 "      Taxon x sample matrix of clade reads (-u: distinct k-mers) to stdout, for taxa with\n"
                 "      >= min_reads in some sample (at -r rank if given). -m: [1]\n"
                 "  index [-a] <out.bidx> <in.bprof>...\n"
                 "      A searchable index of samples: per sample its taxon and clade reads and one k-mer counter\n"
                 "      of all its classified k-mers. -a appends to an existing index.\n"
                 "  search [-n hits] [-M emd|bc|jaccard|containment] <index.bidx> <query.bprof>...\n"
                 "      The indexed samples nearest each query: tree EMD [default], Bray-Curtis (taxa), k-mer\n"
                 "      Jaccard or containment of the query's k-mers in the sample's. -n: [10]\n"
                 "  diff [-r rank] [-m min_reads] [-N permutations] [-s seed] [-K] <groups.tsv>\n"
                 "      Differential abundance between two groups; groups.tsv: \"<profile>\\t<group>\" per line,\n"
                 "      two groups (the first named is A). Per clade (at -r rank, default every clade) present in\n"
                 "      >= 2 samples: Welch t on log proportions, permutation p, Benjamini-Hochberg q; with -K\n"
                 "      also PERMANOVA on the samples' k-mer Jaccard distances for that clade (k-mer content\n"
                 "      shifts, e.g. strain replacement, at equal abundance) and the share of B's pooled k-mers\n"
                 "      absent from A's. Header lines: community PERMANOVA on tree EMD and Bray-Curtis.\n"
                 "      -N: [9999] (all relabelings when fewer), -s: [1], -m: [1]\n", arg);
    return EXIT_FAILURE;
}
uint8_t parse_rank(const char *r) {
    const uint8_t ret = rank_level(r);
    if(ret == RANK_UNKNOWN) RUNTIME_ERROR(std::string("unknown rank ") + r);
    return ret;
}
// Loads the profiles one at a time into summaries, checking they share a database.
template<typename Keep>
std::vector<ProfileSummary> load_summaries(char **paths, int n, Keep &&keep, std::unordered_map<tax_t, uint8_t> &ranks, bool &has_ranks) {
    std::vector<ProfileSummary> ret;
    Profile first;
    has_ranks = true;
    for(int i = 0; i < n; ++i) {
        Profile pr = Profile::read(paths[i]);
        if(i == 0) {
            first.k = pr.k; first.canonical = pr.canonical; first.kind = pr.kind; first.p = pr.p; first.fingerprint = pr.fingerprint;
        } else if(const auto why = first.incompatibility(pr); why.size()) {
            RUNTIME_ERROR(std::string(paths[i]) + " is not comparable with " + paths[0] + ": " + why);
        }
        has_ranks &= bool(pr.has_ranks);
        for(const auto &kv: pr.taxa) if(kv.second.rank != RANK_UNKNOWN) ranks.emplace(kv.first, kv.second.rank);
        ret.push_back(summarize(std::move(pr), keep));
    }
    return ret;
}
// Taxon reads (rank == RANK_UNKNOWN), else clade reads of the taxa at that rank.
std::unordered_map<tax_t, uint64_t> rank_counts(const ProfileSummary &s, uint8_t rank, const std::unordered_map<tax_t, uint8_t> &ranks) {
    if(rank == RANK_UNKNOWN) return s.taxon_reads;
    std::unordered_map<tax_t, uint64_t> ret;
    for(const auto &kv: s.clade_reads) {
        const auto it = ranks.find(kv.first);
        if(it != ranks.end() && it->second == rank && kv.second) ret.emplace(kv.first, kv.second);
    }
    return ret;
}
} // namespace

int profile_main(int argc, char *argv[]) {
    if(argc < 2) return profile_usage(*argv);
    const std::string sub = argv[1];
    --argc, ++argv;
    int c;
    if(sub == "report") {
        std::string dir;
        unsigned b = 16;
        while((c = getopt(argc, argv, "B:b:h?")) >= 0) switch(c) {
            case 'B': dir = optarg; break;
            case 'b': b = unsigned(std::atoi(optarg)); break;
            default: return profile_usage(argv[-1]);
        }
        if(argc - optind < 1 || argc - optind > 2) return profile_usage(argv[-1]);
        Profile pr = Profile::read(argv[optind]);
        if(dir.size() && UniqSketch(pr.kind) != UniqSketch::BBMH) RUNTIME_ERROR("-B needs a b-bit MinHash profile (classify -m or -B)");
        std::unordered_map<tax_t, uint64_t> db;
        for(const auto &kv: pr.taxa) db.emplace(kv.first, kv.second.db_kmers);
        const CompactTaxMap parents = pr.parents();
        auto table = pr.to_table();
        std::FILE *ofp = argc - optind == 2 ? std::fopen(argv[optind + 1], "w"): stdout;
        if(!ofp) RUNTIME_ERROR(std::string("could not open ") + argv[optind + 1]);
        write_uniq_rows(ofp, *table, parents, db, [&](tax_t t) -> std::string {
            return pr.has_ranks ? profile_rank_name(pr.taxa.at(t).rank): "-";
        }, dir.size() ? &dir: nullptr, b);
        if(ofp != stdout) std::fclose(ofp);
        return EXIT_SUCCESS;
    }
    if(sub == "merge") {
        std::string name;
        while((c = getopt(argc, argv, "n:h?")) >= 0) switch(c) {
            case 'n': name = optarg; break;
            default: return profile_usage(argv[-1]);
        }
        if(argc - optind < 2) return profile_usage(argv[-1]);
        const std::string out = argv[optind];
        Profile pr = Profile::read(argv[optind + 1]);
        for(int i = optind + 2; i < argc; ++i) pr.merge(Profile::read(argv[i]));
        if(name.empty()) {
            name = out.substr(out.find_last_of('/') + 1);
            if(endswith(name, ".bprof")) name.resize(name.size() - 6);
        }
        pr.name = name;
        pr.write(out);
        return EXIT_SUCCESS;
    }
    if(sub == "compare" || sub == "table") {
        uint8_t rank = RANK_UNKNOWN;
        uint64_t min_reads = 1;
        std::string kpath;
        bool uniq = false;
        while((c = getopt(argc, argv, "r:m:K:uh?")) >= 0) switch(c) {
            case 'r': rank = parse_rank(optarg); break;
            case 'm': min_reads = std::strtoull(optarg, nullptr, 10); break;
            case 'K': kpath = optarg; break;
            case 'u': uniq = true; break;
            default: return profile_usage(argv[-1]);
        }
        const int n = argc - optind;
        if(n < (sub == "compare" ? 2: 1)) return profile_usage(argv[-1]);
        std::unordered_map<tax_t, uint8_t> ranks;
        bool has_ranks;
        // Only -K keeps clade counters (at the rank, with enough reads).
        auto keep = [&](tax_t, uint64_t reads, uint8_t r) {
            return kpath.size() && reads >= min_reads && (rank == RANK_UNKNOWN || r == rank);
        };
        auto sums = load_summaries(argv + optind, n, keep, ranks, has_ranks);
        if(rank != RANK_UNKNOWN && !has_ranks) RUNTIME_ERROR("-r needs profiles with ranks (classify with a nodes.dmp-style taxonomy)");
        std::vector<std::unordered_map<tax_t, uint64_t>> counts;
        for(const auto &s: sums) counts.push_back(rank_counts(s, rank, ranks));
        auto rank_label = [&](tax_t t) -> std::string {
            if(!has_ranks) return "-";
            const auto it = ranks.find(t);
            return profile_rank_name(it == ranks.end() ? uint8_t(RANK_UNKNOWN): it->second);
        };
        if(sub == "table") {
            std::set<tax_t> rows;
            for(const auto &m: counts) for(const auto &kv: m) if(kv.second >= min_reads) rows.insert(kv.first);
            std::fprintf(stdout, "#taxid\trank");
            for(const auto &s: sums) std::fprintf(stdout, "\t%s", s.name.data());
            std::fputc('\n', stdout);
            for(const tax_t t: rows) {
                std::fprintf(stdout, "%u\t%s", t, rank_label(t).data());
                for(size_t i = 0; i < sums.size(); ++i) {
                    if(uniq) {
                        const auto it = sums[i].clade_uniq.find(t);
                        std::fprintf(stdout, "\t%.0f", it == sums[i].clade_uniq.end() ? 0.: it->second);
                    } else {
                        const auto it = counts[i].find(t);
                        std::fprintf(stdout, "\t%llu", it == counts[i].end() ? 0ull: (unsigned long long)it->second);
                    }
                }
                std::fputc('\n', stdout);
            }
            return EXIT_SUCCESS;
        }
        std::fprintf(stdout, "#sample_a\tsample_b\tclassified_a\tclassified_b\tbray_curtis\tjaccard_dist\tshared_taxa\ttree_emd\n");
        for(size_t i = 0; i < sums.size(); ++i) for(size_t j = i + 1; j < sums.size(); ++j) {
            const auto d = abundance_distance(counts[i], counts[j], sums[i].clade_reads, sums[i].classified,
                                              sums[j].clade_reads, sums[j].classified, min_reads);
            std::fprintf(stdout, "%s\t%s\t%llu\t%llu\t%.6f\t%.6f\t%zu\t%.6f\n", sums[i].name.data(), sums[j].name.data(),
                         (unsigned long long)sums[i].classified, (unsigned long long)sums[j].classified,
                         d.bray_curtis, d.jaccard, d.shared, d.tree_emd);
        }
        if(kpath.size()) {
            std::FILE *kfp = std::fopen(kpath.data(), "w");
            if(!kfp) RUNTIME_ERROR("could not open " + kpath);
            std::fprintf(kfp, "#taxid\trank\tsample_a\tsample_b\treads_a\treads_b\tuniq_a\tuniq_b\tshared\tjaccard\tcont_a_in_b\tcont_b_in_a\n");
            for(size_t i = 0; i < sums.size(); ++i) for(size_t j = i + 1; j < sums.size(); ++j) {
                std::vector<tax_t> both;
                for(const auto &kv: sums[i].counters) if(sums[j].counters.count(kv.first)) both.push_back(kv.first);
                std::sort(both.begin(), both.end());
                for(const tax_t t: both) {
                    const auto o = counter_overlap(sums[i].counters.at(t), sums[j].counters.at(t));
                    std::fprintf(kfp, "%u\t%s\t%s\t%s\t%llu\t%llu\t%.0f\t%.0f\t%.0f\t%.6f\t%.6f\t%.6f\n", t, rank_label(t).data(),
                                 sums[i].name.data(), sums[j].name.data(),
                                 (unsigned long long)sums[i].clade_reads.at(t), (unsigned long long)sums[j].clade_reads.at(t),
                                 o.a, o.b, o.shared, o.jaccard, o.a > 0. ? std::min(1., o.shared / o.a): 0., o.b > 0. ? std::min(1., o.shared / o.b): 0.);
                }
            }
            std::fclose(kfp);
        }
        return EXIT_SUCCESS;
    }
    if(sub == "index") {
        bool append = false;
        while((c = getopt(argc, argv, "ah?")) >= 0) switch(c) {
            case 'a': append = true; break;
            default: return profile_usage(argv[-1]);
        }
        if(argc - optind < 2) return profile_usage(argv[-1]);
        const std::string out = argv[optind];
        Profile header;
        bool have = append && isfile(out);
        if(have) header = std::move(SampleIndex::read(out, false).header);
        gzFile fp = Profile::open(out, have ? "ab": "wb");
        for(int i = optind + 1; i < argc; ++i) {
            Profile pr = Profile::read(argv[i]);
            if(!have) {
                header.k = pr.k; header.canonical = pr.canonical; header.kind = pr.kind; header.p = pr.p; header.fingerprint = pr.fingerprint;
                SampleIndex::write_header(fp, header);
                have = true;
            } else if(const auto why = header.incompatibility(pr); why.size()) {
                gzclose(fp);
                RUNTIME_ERROR(std::string(argv[i]) + " does not match the index: " + why);
            }
            IndexEntry e = index_entry(summarize(std::move(pr), [](tax_t, uint64_t, uint8_t) {return false;}));
            SampleIndex::write_entry(fp, e);
        }
        if(gzclose(fp) != Z_OK) RUNTIME_ERROR("could not finish index " + out);
        return EXIT_SUCCESS;
    }
    if(sub == "search") {
        size_t nhits = 10;
        std::string metric = "emd";
        while((c = getopt(argc, argv, "n:M:h?")) >= 0) switch(c) {
            case 'n': nhits = std::strtoull(optarg, nullptr, 10); break;
            case 'M': metric = optarg; break;
            default: return profile_usage(argv[-1]);
        }
        if(metric != "emd" && metric != "bc" && metric != "jaccard" && metric != "containment") RUNTIME_ERROR("-M: emd, bc, jaccard or containment");
        if(argc - optind < 2) return profile_usage(argv[-1]);
        SampleIndex idx = SampleIndex::read(argv[optind]);
        std::fprintf(stdout, "#query\thit\tsample\tbray_curtis\ttree_emd\tkmer_jaccard\tcont_query_in_sample\tcont_sample_in_query\n");
        for(int qi = optind + 1; qi < argc; ++qi) {
            Profile pr = Profile::read(argv[qi]);
            if(const auto why = idx.header.incompatibility(pr); why.size()) RUNTIME_ERROR(std::string(argv[qi]) + " does not match the index: " + why);
            ProfileSummary q = summarize(std::move(pr), [](tax_t, uint64_t, uint8_t) {return false;});
            struct Hit {size_t i; AbundanceDistance d; CounterOverlap o; double key;};
            std::vector<Hit> hits;
            for(size_t i = 0; i < idx.entries.size(); ++i) {
                IndexEntry &e = idx.entries[i];
                Hit h{i, abundance_distance(q.taxon_reads, e.taxon_reads, q.clade_reads, q.classified, e.clade_reads, e.classified, 1),
                      counter_overlap(q.all, e.sample), 0.};
                h.key = metric == "emd" ? h.d.tree_emd: metric == "bc" ? h.d.bray_curtis
                      : metric == "jaccard" ? -h.o.jaccard: -(h.o.a > 0. ? h.o.shared / h.o.a: 0.);
                hits.push_back(h);
            }
            std::sort(hits.begin(), hits.end(), [&](const Hit &a, const Hit &b) {return a.key != b.key ? a.key < b.key: idx.entries[a.i].name < idx.entries[b.i].name;});
            for(size_t r = 0; r < std::min(nhits, hits.size()); ++r) {
                const Hit &h = hits[r];
                std::fprintf(stdout, "%s\t%zu\t%s\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n", q.name.data(), r + 1, idx.entries[h.i].name.data(),
                             h.d.bray_curtis, h.d.tree_emd, h.o.jaccard, h.o.a > 0. ? std::min(1., h.o.shared / h.o.a): 0.,
                             h.o.b > 0. ? std::min(1., h.o.shared / h.o.b): 0.);
            }
        }
        return EXIT_SUCCESS;
    }
    if(sub == "diff") {
        uint8_t rank = RANK_UNKNOWN;
        uint64_t min_reads = 1, seed = 1;
        size_t nperm = 9999;
        bool kmers = false;
        while((c = getopt(argc, argv, "r:m:N:s:Kh?")) >= 0) switch(c) {
            case 'r': rank = parse_rank(optarg); break;
            case 'm': min_reads = std::strtoull(optarg, nullptr, 10); break;
            case 'N': nperm = std::strtoull(optarg, nullptr, 10); break;
            case 's': seed = std::strtoull(optarg, nullptr, 10); break;
            case 'K': kmers = true; break;
            default: return profile_usage(argv[-1]);
        }
        if(argc - optind != 1) return profile_usage(argv[-1]);
        std::vector<std::string> paths, gnames;
        std::vector<uint8_t> labels;
        {
            std::ifstream is(argv[optind]);
            if(!is) RUNTIME_ERROR(std::string("could not open ") + argv[optind]);
            for(std::string line; std::getline(is, line);) {
                if(line.empty() || line[0] == '#') continue;
                const size_t tab = line.find('\t');
                if(tab == std::string::npos) RUNTIME_ERROR("groups file: \"<profile>\\t<group>\" per line: " + line);
                const std::string g = line.substr(tab + 1);
                auto it = std::find(gnames.begin(), gnames.end(), g);
                if(it == gnames.end()) {gnames.push_back(g); it = gnames.end() - 1;}
                paths.push_back(line.substr(0, tab));
                labels.push_back(uint8_t(it - gnames.begin()));
            }
        }
        if(gnames.size() != 2) RUNTIME_ERROR("diff compares exactly two groups; the groups file names " + std::to_string(gnames.size()));
        const size_t n = paths.size();
        std::vector<char *> cpaths;
        for(auto &p: paths) cpaths.push_back(p.data());
        std::unordered_map<tax_t, uint8_t> ranks;
        bool has_ranks;
        auto keep = [&](tax_t, uint64_t reads, uint8_t r) {return kmers && reads >= min_reads && (rank == RANK_UNKNOWN || r == rank);};
        auto sums = load_summaries(cpaths.data(), int(n), keep, ranks, has_ranks);
        if(rank != RANK_UNKNOWN && !has_ranks) RUNTIME_ERROR("-r needs profiles with ranks (classify with a nodes.dmp-style taxonomy)");
        auto rank_label = [&](tax_t t) -> std::string {
            if(!has_ranks) return "-";
            const auto it = ranks.find(t);
            return profile_rank_name(it == ranks.end() ? uint8_t(RANK_UNKNOWN): it->second);
        };
        size_t na = std::count(labels.begin(), labels.end(), 0);
        std::fprintf(stdout, "#groups\tA=%s\t%zu\tB=%s\t%zu\tseed=%llu\tpermutations=%zu\n", gnames[0].data(), na, gnames[1].data(), n - na,
                     (unsigned long long)seed, nperm);
        // community-wide
        std::vector<std::unordered_map<tax_t, uint64_t>> counts;
        for(const auto &s: sums) counts.push_back(rank_counts(s, rank, ranks));
        std::vector<double> demd(n * n, 0.), dbc(n * n, 0.);
        for(size_t i = 0; i < n; ++i) for(size_t j = i + 1; j < n; ++j) {
            const auto d = abundance_distance(counts[i], counts[j], sums[i].clade_reads, sums[i].classified, sums[j].clade_reads, sums[j].classified, min_reads);
            demd[i * n + j] = demd[j * n + i] = d.tree_emd;
            dbc[i * n + j] = dbc[j * n + i] = d.bray_curtis;
        }
        for(const auto &m: {std::make_pair("tree_emd", &demd), std::make_pair("bray_curtis", &dbc)}) {
            const double f = permanova_f(*m.second, n, labels);
            const auto pp = permutation_p(labels, [&](const std::vector<uint8_t> &l) {return permanova_f(*m.second, n, l);}, nperm, seed);
            std::fprintf(stdout, "#permanova\t%s\tF=%.6g\tp=%.6g\trelabelings=%zu\n", m.first, f, pp.first, pp.second);
        }
        // per clade
        std::set<tax_t> taxa;
        for(const auto &m: counts) for(const auto &kv: m) taxa.insert(kv.first);
        if(rank == RANK_UNKNOWN) for(const auto &s: sums) for(const auto &kv: s.clade_reads) taxa.insert(kv.first);
        struct Row {tax_t t; size_t prev[2]; double mean[2], lfc, p, kp = NAN, kf = NAN, within = NAN, between = NAN, novel = NAN; size_t krel = 0;};
        std::vector<Row> rows;
        for(const tax_t t: taxa) {
            Row r{t, {0, 0}, {0., 0.}, 0., NAN};
            std::vector<double> x(n);
            for(size_t i = 0; i < n; ++i) {
                const auto it = sums[i].clade_reads.find(t);
                const uint64_t v = it == sums[i].clade_reads.end() ? 0: it->second;
                r.prev[labels[i]] += v >= min_reads;
                r.mean[labels[i]] += sums[i].classified ? double(v) / sums[i].classified: 0.;
                x[i] = std::log((v + .5) / (sums[i].classified + .5));   // 0 iff the clade holds every read
            }
            if(r.prev[0] + r.prev[1] < 2) continue;
            r.mean[0] /= na; r.mean[1] /= (n - na);
            double mx[2] = {0., 0.};
            for(size_t i = 0; i < n; ++i) mx[labels[i]] += x[i];
            r.lfc = (mx[1] / (n - na) - mx[0] / na) / std::log(2.);
            r.p = permutation_p(labels, [&](const std::vector<uint8_t> &l) {return welch_abs_t(x, l);}, nperm, seed ^ uniq_hash(uint64_t(t))).first;
            if(kmers) {
                std::vector<size_t> who;
                for(size_t i = 0; i < n; ++i) if(sums[i].counters.count(t)) who.push_back(i);
                std::vector<uint8_t> kl;
                for(const size_t i: who) kl.push_back(labels[i]);
                const size_t kb = std::count(kl.begin(), kl.end(), 1);
                if(kb >= 1 && who.size() - kb >= 1) {
                    const size_t m = who.size();
                    std::vector<double> dj(m * m, 0.);
                    double ws = 0., bs = 0.; size_t wn = 0, bn = 0;
                    for(size_t a = 0; a < m; ++a) for(size_t b = a + 1; b < m; ++b) {
                        const double d = 1. - counter_overlap(sums[who[a]].counters.at(t), sums[who[b]].counters.at(t)).jaccard;
                        dj[a * m + b] = dj[b * m + a] = d;
                        if(kl[a] == kl[b]) ws += d, ++wn; else bs += d, ++bn;
                    }
                    r.within = wn ? 1. - ws / wn: NAN; r.between = bn ? 1. - bs / bn: NAN;
                    if(kb >= 2 && m - kb >= 2) {
                        r.kf = permanova_f(dj, m, kl);
                        const auto pp = permutation_p(kl, [&](const std::vector<uint8_t> &l) {return permanova_f(dj, m, l);}, nperm, seed ^ uniq_hash(uint64_t(t) + 1));
                        r.kp = pp.first; r.krel = pp.second;
                    }
                    UniqCounter pool[2];
                    bool init[2] = {false, false};
                    for(size_t a = 0; a < m; ++a) {
                        UniqCounter &c = sums[who[a]].counters.at(t);
                        if(!init[kl[a]]) pool[kl[a]] = c.clone(), init[kl[a]] = true;
                        else pool[kl[a]].merge(c);
                    }
                    const auto o = counter_overlap(pool[0], pool[1]);
                    r.novel = o.b > 0. ? std::max(0., 1. - o.shared / o.b): NAN;
                }
            }
            rows.push_back(r);
        }
        std::vector<double> ps, kps;
        for(const auto &r: rows) ps.push_back(r.p), kps.push_back(r.kp);
        const auto qs = bh_qvalues(ps), kqs = bh_qvalues(kps);
        std::vector<size_t> order(rows.size());
        std::iota(order.begin(), order.end(), size_t(0));
        std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
            const double pa = std::isnan(rows[a].p) ? 2.: rows[a].p, pb = std::isnan(rows[b].p) ? 2.: rows[b].p;
            return pa != pb ? pa < pb: rows[a].t < rows[b].t;
        });
        std::fprintf(stdout, "#taxid\trank\tprev_a\tprev_b\tmean_prop_a\tmean_prop_b\tlog2fc\tp\tq");
        if(kmers) std::fprintf(stdout, "\tkmer_jaccard_within\tkmer_jaccard_between\tkmer_F\tkmer_p\tkmer_q\tkmer_relabelings\tb_kmers_novel_vs_a");
        std::fputc('\n', stdout);
        for(const size_t i: order) {
            const Row &r = rows[i];
            std::fprintf(stdout, "%u\t%s\t%zu\t%zu\t%.6g\t%.6g\t%.4f\t%.6g\t%.6g", r.t, rank_label(r.t).data(), r.prev[0], r.prev[1],
                         r.mean[0], r.mean[1], r.lfc, r.p, qs[i]);
            if(kmers) std::fprintf(stdout, "\t%.6f\t%.6f\t%.6g\t%.6g\t%.6g\t%zu\t%.6f", r.within, r.between, r.kf, r.kp, kqs[i], r.krel, r.novel);
            std::fputc('\n', stdout);
        }
        return EXIT_SUCCESS;
    }
    return profile_usage(argv[-1]);
}

int main(int argc, char *argv[]) {
    const static std::unordered_map<std::string, bns::MainFnPtr> md {
        {"phase1",   phase1_main},
        {"prebuild", phase1_main},
        {"p1",       phase1_main},
        {"phase2",   phase2_main},
        {"build",    phase2_main},
        {"p2",       phase2_main},
        {"lca",      phase1_main},
        {"hist",     hist_main},
        {"metatree", metatree_main},
        {"export",   export_main},
        {"classify", classify_main},
        {"compact",  compact_main},
        {"profile",  profile_main}
    };
    if(std::find_if(argv, argv + argc, [&](char *s) {return std::strcmp("-v", s) == 0 || std::strcmp("--version", s) == 0;}) != argv + argc) {
        std::fprintf(stdout, "bonsai|%s\n", BONSAI_VERSION);
        return EXIT_SUCCESS;
    }
    // Report errors (RUNTIME_ERROR and friends) and exit 1 instead of aborting.
    try {
        return (argc > 1 && md.find(argv[1]) != md.end() ? md.find(argv[1])->second
                                                         : err_main)(argc - 1, argv + 1);
    } catch(const std::exception &e) {
        std::fprintf(stderr, "bonsai: error: %s\n", e.what());
        return EXIT_FAILURE;
    }
}
