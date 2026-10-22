#ifndef __BNS_KRAKEN_EXPORT_H
#define __BNS_KRAKEN_EXPORT_H

/*
Export a (possibly grafted) PlacementTree as taxonomy files that
Kraken 2 and KrakenUniq can ingest directly, so downstream users can
build classifiers on the refined tree without running bonsai's own.

Layout produced under the output directory:
  nodes.dmp        NCBI nodes.dmp-style: "<id>\t|\t<parent>\t|\t<rank>\t|"
  names.dmp        "<id>\t|\t<name>\t|\t\t|\tscientific name\t|"
  seqid2taxid.map          "<seqid>\t<taxid>" per genome, keyed on the
                           ORIGINAL header tokens (KrakenUniq input)
  seqid2taxid_kraken2.map  same, keyed on the rewritten tokens
                           "seqid|kraken:taxid|taxid" -- copy to
                           "<k2db>/seqid2taxid.map" (Kraken 2's build_db
                           looks the full rewritten header token up in it)
  library/         the input fastas, decompressed, as <stem>.fna
                   (KrakenUniq's build finds *.fna/*.fa/... and concatenates)
  library_kraken2/ <stem>.fna files with headers rewritten to
                   "><seqid>|kraken:taxid|<taxid> ..." (Kraken 2 custom
                   library; copy the CONTENTS into the db library dir)
Kraken 2's parser (src/taxonomy.cc) splits on "\t|\t", reads node id,
parent id, and rank from the first three fields of nodes.dmp, records a
name only when the names.dmp entry's class is "scientific name", and
forces node 1's parent to 0. KrakenUniq (scripts/build_db.sh) expects
library fastas, taxonomy/nodes.dmp, taxonomy/names.dmp, and
seqid2taxid.map.

After exporting:
  # KrakenUniq
  mkdir -p KUDB && cp -r <outdir>/library <outdir>/seqid2taxid.map KUDB/
  mkdir -p KUDB/taxonomy && cp <outdir>/nodes.dmp <outdir>/names.dmp KUDB/taxonomy/
  krakenuniq-build --db KUDB
  # Kraken 2
  mkdir -p K2DB/taxonomy K2DB/library
  cp <outdir>/nodes.dmp <outdir>/names.dmp K2DB/taxonomy/
  cp <outdir>/library_kraken2/. K2DB/library/
  kraken2-build --build --db K2DB
*/

#include "hashmap.h"
#include "placement.h"

#include <sys/stat.h>

namespace bns {

struct KrakenExportResult {
    std::vector<std::pair<std::string, tax_t>> seqid2taxid; // seqid -> exported leaf taxid
    tax_t max_taxid;                                        // highest id used
    u32 ngrafted;                                           // synthetic nodes exported
};

class KrakenExporter {
public:
    static KrakenExportResult write(const PlacementTree &tree,
                     const std::vector<std::string> &seqids,       // one per genome, in tree order
                     const std::vector<std::string> &genome_paths, // fasta for each genome
                     const char *outdir,
                     const char *graft_prefix = "grafted_",
                     const bns::umap<tax_t, std::string> *taxon_names = nullptr,
                     const bns::umap<tax_t, std::string> *taxon_ranks = nullptr,
                     std::FILE *log = stderr) {
        // Assign export ids: taxonomy nodes keep their ids, leaves and
        // synthetic nodes get fresh ids above every taxonomy id.
        tax_t max_taxid(1);
        for(u32 i(0); i < tree.nnodes(); ++i)
            if(tree.is_alive(i) && tree.node_taxid(i) > max_taxid) max_taxid = tree.node_taxid(i);
        std::vector<tax_t> xid(tree.nnodes(), 0);
        tax_t next(max_taxid + 1);
        std::function<void(u32)> assign = [&](u32 x) {
            if(!tree.is_alive(x)) return;
            if(tree.node_taxid(x)) xid[x] = tree.node_taxid(x);
            else xid[x] = next++;
            for(const auto c: tree.children_of(x)) assign(c);
        };
        assign(tree.root());
        max_taxid = next - 1;

        // Names and ranks.
        auto name_of = [&](u32 x) -> std::string {
            if(x == tree.root()) return "root";
            if(tree.node_taxid(x)) {
                if(taxon_names) {
                    const auto it(taxon_names->find(tree.node_taxid(x)));
                    if(it != taxon_names->end()) return it->second;
                }
                return "taxon_" + std::to_string(tree.node_taxid(x));
            }
            if(tree.is_synthetic(x)) return std::string(graft_prefix) + "node_" + std::to_string(xid[x]);
            for(u32 i(0); i < seqids.size(); ++i) if(tree.leaf_node(i) == x) return "genome_" + seqids[i];
            return "node_" + std::to_string(xid[x]);
        };
        auto rank_of = [&](u32 x) -> std::string {
            if(tree.node_taxid(x) && taxon_ranks) {
                const auto it(taxon_ranks->find(tree.node_taxid(x)));
                if(it != taxon_ranks->end()) return it->second;
            }
            return "no rank";
        };

        // Write nodes.dmp / names.dmp.
        const std::string npath(std::string(outdir) + "/nodes.dmp"), mpath(std::string(outdir) + "/names.dmp");
        std::FILE *nfp(std::fopen(npath.data(), "w")), *mfp(std::fopen(mpath.data(), "w"));
        if(!nfp || !mfp) RUNTIME_ERROR(ks::sprintf("Could not open taxonomy files for writing in %s", outdir).data());
        std::function<void(u32)> emit = [&](u32 x) {
            if(!tree.is_alive(x)) return;
            const tax_t parent(x == tree.root() ? xid[x] : xid[tree.parent_of(x)]);
            std::fprintf(nfp, "%u\t|\t%u\t|\t%s\t|\n", xid[x], parent, rank_of(x).data());
            std::fprintf(mfp, "%u\t|\t%s\t|\t\t|\tscientific name\t|\n", xid[x], name_of(x).data());
            for(const auto c: tree.children_of(x)) emit(c);
        };
        emit(tree.root());
        std::fclose(nfp);
        std::fclose(mfp);

        // seqid2taxid.map and libraries. KrakenUniq's set_lcas keys on the
        // original header tokens (its library is verbatim), while kraken2's
        // build_db keys on the full rewritten header token
        // ("seqid|kraken:taxid|taxid"), so each library gets its own map.
        KrakenExportResult ret;
        ret.max_taxid = max_taxid;
        ret.ngrafted = tree.nsynthetic();
        std::FILE *sfp(std::fopen((std::string(outdir) + "/seqid2taxid.map").data(), "w"));
        std::FILE *s2fp(std::fopen((std::string(outdir) + "/seqid2taxid_kraken2.map").data(), "w"));
        if(!sfp || !s2fp) RUNTIME_ERROR("Could not open seqid2taxid map for writing");
        const std::string ldir(std::string(outdir) + "/library"), l2dir(std::string(outdir) + "/library_kraken2");
        ::mkdir(ldir.data(), 0755);
        ::mkdir(l2dir.data(), 0755);
        // Stems from different directories can collide (a/x.fa, b/x.fa);
        // dedupe so one genome's fasta never silently overwrites another's.
        bns::uset<std::string> used_stems;
        for(u32 i(0); i < seqids.size(); ++i) {
            const tax_t leaf(xid[tree.leaf_node(i)]);
            std::fprintf(sfp, "%s\t%u\n", seqids[i].data(), leaf);
            std::fprintf(s2fp, "%s|kraken:taxid|%u\t%u\n", seqids[i].data(), leaf, leaf);
            ret.seqid2taxid.emplace_back(seqids[i], leaf);
            // KrakenUniq's build finds library files by extension
            // (*.fna/*.fa/...; no .gz) and concatenates them verbatim, and
            // kraken2's build_kraken2_db.sh matches only '*.fna'/'*.faa';
            // both read the copies decompressed, so emit '.fna' stems.
            const std::string base(library_stem(genome_paths[i]));
            std::string stem(base);
            for(u32 n(2); !used_stems.insert(stem).second; stem = base + "_" + std::to_string(n++));
            copy_fasta(genome_paths[i], (ldir + "/" + stem + ".fna").data(), nullptr, 0);
            copy_fasta(genome_paths[i], (l2dir + "/" + stem + ".fna").data(), seqids[i].data(), leaf);
        }
        std::fclose(sfp);
        std::fclose(s2fp);
        if(log) {
            std::fprintf(log, "Exported %zu genomes (%zu taxonomy nodes kept, %u grafted nodes, max taxid %u) to %s\n",
                         seqids.size(), size_t(tree.nnodes()), ret.ngrafted, ret.max_taxid, outdir);
            std::fprintf(log, "KrakenUniq:  mkdir -p KUDB/taxonomy && cp -r %s/library %s/seqid2taxid.map KUDB/ && cp %s/nodes.dmp %s/names.dmp KUDB/taxonomy/ && krakenuniq-build --db KUDB\n",
                         outdir, outdir, outdir, outdir);
            std::fprintf(log, "Kraken 2:    mkdir -p K2DB/taxonomy K2DB/library && cp %s/nodes.dmp %s/names.dmp K2DB/taxonomy/ && cp %s/seqid2taxid_kraken2.map K2DB/seqid2taxid.map && cp %s/library_kraken2/* K2DB/library/ && kraken2-build --build --db K2DB\n",
                         outdir, outdir, outdir, outdir);
        }
        return ret;
    }

private:
    static std::string base_name(const std::string &path) {
        const auto pos(path.find_last_of('/'));
        return pos == std::string::npos ? path : path.substr(pos + 1);
    }
    // Library file stem: basename with any compression suffix stripped and
    // the '.fna'/'.fa'/... extension removed; kraken2's build_kraken2_db.sh
    // matches only '*.fna'/'*.faa' in the library and krakenuniq's find has
    // no '.gz', so the caller appends a fresh '.fna'.
    static std::string library_stem(const std::string &path) {
        std::string ret(base_name(path));
        if(ret.size() > 3 && ret.compare(ret.size() - 3, 3, ".gz") == 0) ret.resize(ret.size() - 3);
        else if(ret.size() > 5 && ret.compare(ret.size() - 5, 5, ".gzip") == 0) ret.resize(ret.size() - 5);
        static const char *exts[] {".fna", ".faa", ".fa", ".ffn", ".fsa", ".fasta"};
        for(const auto ext: exts) {
            const size_t el(std::strlen(ext));
            if(ret.size() > el && ret.compare(ret.size() - el, el, ext) == 0) {
                ret.resize(ret.size() - el);
                break;
            }
        }
        return ret;
    }
    static ssize_t gz_read_line(gzFile fp, char **buf, size_t *cap) {
        if(*cap == 0) {*cap = 4096; *buf = static_cast<char *>(std::malloc(*cap));}
        size_t n(0);
        int c(-1);
        while((c = gzgetc(fp)) != -1) {
            if(n + 1 >= *cap) {*cap <<= 1; *buf = static_cast<char *>(std::realloc(*buf, *cap));}
            (*buf)[n++] = c;
            if(c == '\n') break;
        }
        if(c == -1 && n == 0) return -1;
        (*buf)[n] = 0;
        return n;
    }
    // Copy a fasta, optionally rewriting each header's first token to
    // "<seqid>|kraken:taxid|<taxid>".
    static void copy_fasta(const std::string &src, const char *dst, const char *seqid, tax_t taxid) {
        gzFile in(gzopen(src.data(), "rb"));
        if(!in) RUNTIME_ERROR(ks::sprintf("Could not open %s", src.data()).data());
        std::FILE *out(std::fopen(dst, "w"));
        if(!out) RUNTIME_ERROR(ks::sprintf("Could not open %s", dst).data());
        char *line(nullptr);
        size_t cap(0);
        while(gz_read_line(in, &line, &cap) >= 0) {
            if(line[0] == '>' && seqid) {
                char *rest(line + 1);
                while(*rest == ' ') ++rest;
                while(*rest && !std::isspace(*rest)) ++rest; // skip original first token
                std::fprintf(out, ">%s|kraken:taxid|%u%s", seqid, taxid, rest);
            } else std::fputs(line, out);
        }
        std::free(line);
        gzclose(in);
        std::fclose(out);
    }
};

} // namespace bns

#endif // __BNS_KRAKEN_EXPORT_H
