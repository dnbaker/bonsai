#include "util.h"
#include <set>
#include <cassert>
#include <cstring>
#include <zlib.h>
#include <sstream>
#include <fstream>

namespace emp {

std::size_t count_lines(const char *fn) noexcept {
    std::FILE *fp(std::fopen(fn, "r"));
    std::size_t bufsz = 4096;
    char *buf((char *)std::malloc(bufsz));
    ssize_t len;
    std::size_t n(0);
    while((len = getline(&buf, &bufsz, fp)) >= 0) ++n;
    std::free(buf);
    std::fclose(fp);
    return n;
}
std::unordered_map<tax_t, std::vector<tax_t>> invert_parent_map(khash_t(p) *tax) noexcept {
    std::unordered_map<tax_t, std::vector<tax_t>> ret;
    typename std::unordered_map<tax_t, std::vector<tax_t>>::iterator m;
    for(khiter_t ki(0); ki < kh_end(tax); ++ki) {
        if(!kh_exist(tax, ki)) continue;
        if((m = ret.find(kh_val(tax, ki))) == ret.end()) m = ret.emplace(kh_val(tax, ki), std::vector<tax_t>{kh_key(tax, ki)}).first;
        else                                             m->second.emplace_back(kh_key(tax, ki));
    }
#if !NDEBUG
    for(auto &pair: ret) {
        for(auto child: pair.second) assert(get_parent(tax, child) == pair.first);
    }
#endif
    return ret;
}

std::vector<tax_t> get_all_descendents(const std::unordered_map<tax_t, std::vector<tax_t>> &map, tax_t tax) {
    std::vector<tax_t> ret;
    auto m(map.find(tax));
    if(m != map.end()) {
        LOG_DEBUG("Found a match for parent %u with %zu children.\n", tax, m->second.size());
        for(tax_t newtax: m->second) {
            assert(newtax != tax);
            ret.push_back(newtax);
            auto desc(get_all_descendents(map, newtax));
            for(auto d: desc) std::fprintf(stderr, "Descendent %u of total %zu descendences of %u\n", d, desc.size(), tax);
            ret.insert(ret.end(), desc.begin(), desc.end());
        }
    }
    LOG_DEBUG("Returning descendents. Size: %zu\n", ret.size());
    return ret;
}

// Rewritten from Kraken's source code.
// Consider rewriting to use kbtree instead of std::map.
// lh3's benchmarks indicate that his is only as fast as std::map,
// though it uses less than half the memory. The thing is, taxonomic trees
// aren't that deep.
// I don't think we'd gain anything practical with that.
tax_t lca(const khash_t(p) *map, tax_t a, tax_t b) noexcept {
    // Use std::set to use RB trees for small set rather than hash table.
    std::set<tax_t> nodes;
    if(a == b) return a;
    if(b == 0) return a;
    if(a == 0) return b;
    khint_t ki;
    while(a) {
        nodes.insert(a);
        if((ki = kh_get(p, map, a)) == kh_end(map)) {
            fprintf(stderr, "Missing taxid %u. Returning -1!\n", a);
            return (tax_t)-1;
        }
        a = kh_val(map, ki);
    }
    while(b) {
        if(nodes.find(b) != nodes.end()) return b;
        if((ki = kh_get(p, map, b)) == kh_end(map)) {
            fprintf(stderr, "Missing taxid %u. Returning -1!\n", b);
            return (tax_t)-1;
        }
        b = kh_val(map, ki);
    }
    return 1;
}

unsigned node_dist(khash_t(p) *map, tax_t leaf, tax_t root) noexcept {
    unsigned ret(0);
    khint_t ki;
    while(leaf) {
        if((ki = kh_get(p, map, leaf)) == kh_end(map)) {
            fprintf(stderr, "Tax ID %u missing. Abort!\n", leaf);
            std::exit(1);
        }
        ++ret;
        if((leaf = kh_val(map, ki)) == root) return ret;
    }
    LOG_EXIT("leaf %u is not a child of root %u\n", leaf, root);
    return ret;
}
unsigned node_depth(khash_t(p) *map, tax_t a) noexcept {
    unsigned ret(0);
    khint_t ki;
    while(a) {
        if((ki = kh_get(p, map, a)) == kh_end(map)) {
            fprintf(stderr, "Tax ID %u missing. Abort!\n", a);
            std::exit(1);
        }
        a = kh_val(map, ki);
        ++ret;
    }
    return ret;
}

khash_t(name) *build_name_hash(const char *fn) noexcept {
    LOG_INFO("Parsing name hash from %s\n", fn);
    std::size_t bufsz(2048), namelen;
    char *buf((char *)std::malloc(bufsz));
    ssize_t len;
    std::FILE *fp(std::fopen(fn, "r"));
    khash_t(name) *ret(kh_init(name));
    kh_resize(name, ret, count_lines(fn));
    char *p;
    int khr;
    khint_t ki;
    while((len = getline(&buf, &bufsz, fp)) >= 0) {
        switch(*buf) case '\0': case '\n': case '#': continue;
        p = strchr(buf, '\t');
        *p = 0;
        ki = kh_put(name, ret, buf, &khr);
        if(khr == 0) { // Key already present.
            LOG_INFO("Key %s already present. Updating value from %i to %s\n", kh_key(ret, ki), kh_val(ret, ki), p + 1);
        }  else {
            namelen = p - buf;
            kh_key(ret, ki) = (char *)std::malloc(namelen + 1);
            memcpy((void*)kh_key(ret, ki), buf, namelen);
           ((char *)kh_key(ret, ki))[namelen] = '\0';
        }
        kh_val(ret, ki) = atoi(p + 1);
    }
    std::fclose(fp);
    std::free(buf);
    return ret;
}

void destroy_name_hash(khash_t(name) *hash) noexcept {
    for(khint_t ki(kh_begin(hash)); ki != kh_end(hash); ++ki)
        if(kh_exist(hash, ki))
            std::free((void *)kh_key(hash, ki));
    kh_destroy(name, hash);
}

std::map<tax_t, tax_t> build_kraken_tax(const std::string &fname) {
    const char *fn(fname.data());
    std::FILE *fp(std::fopen(fn, "r"));
    std::map<tax_t, tax_t> ret;
    std::size_t bufsz = 4096;
    char *buf((char *)std::malloc(bufsz));
    ssize_t len;
    tax_t t1, t2;
    while((len = getline(&buf, &bufsz, fp)) >= 0) {
        switch(*buf) case '\n': case '0': case '#': continue;
        t1 = atoi(buf);
        t2 = atoi(strchr(buf, '|') + 2);
        ret[t1] = t2;
    }
    ret[1] = 0;
    std::fclose(fp);
    std::free(buf);
    return ret;
}
uint32_t lca(std::map<uint32_t, uint32_t> &parent_map, uint32_t a, uint32_t b)
{
  if (a == 0 || b == 0)
    return a ? a : b;

  std::set<uint32_t> a_path;
  while (a) {
    a_path.insert(a);
    a = parent_map[a];
  }
  while (b) {
    if (a_path.count(b) > 0)
      return b;
    b = parent_map[b];
  }
  return 1;
}

khash_t(p) *build_parent_map(const char *fn) noexcept {
    std::FILE *fp(std::fopen(fn, "r"));
    khash_t(p) *ret(kh_init(p));
    std::size_t bufsz = 4096;
    char *buf((char *)std::malloc(bufsz));
    ssize_t len;
    khint_t ki;
    int khr;
    while((len = getline(&buf, &bufsz, fp)) >= 0) {
        switch(*buf) case '\n': case '0': case '#': continue;
        ki = kh_put(p, ret, atoi(buf), &khr);
        kh_val(ret, ki) = atoi(strchr(buf, '|') + 2);
    }
    ki = kh_put(p, ret, 1, &khr);
    kh_val(ret, ki) = 0; // Root of the tree.
    std::fclose(fp);
    std::free(buf);
    return ret;
}

void kset_union(khash_t(all) *a, khash_t(all) *b) noexcept {
    khint_t ki2;
    int khr;
    for(ki2 = kh_begin(b); ki2 != kh_end(b); ++ki2)
        if(kh_exist(b, ki2))
            kh_put(all, a, kh_key(b, ki2), &khr);
}

 // Tree resolution: take all hit taxa (plus ancestors), then
 // return leaf of highest weighted leaf-to-root path.
 // Taken, slightly modified from Kraken
tax_t resolve_tree(std::map<tax_t, tax_t> &hit_counts,
                      khash_t(p) *parent_map) noexcept
{
  std::set<tax_t> max_taxa;
  tax_t max_taxon(0), max_score(0);

  // Sum each taxon's LTR path
  for(auto it(hit_counts.cbegin()), e(hit_counts.cend()); it != e; ++it) {
    tax_t taxon(it->first), node(taxon), score(0);
    khiter_t ki;
    // Instead of while node > 0
    while(node) {
        score += hit_counts[node];
        ki = kh_get(p, parent_map, node);
        node = kh_val(parent_map, ki);
    }
    if (score > max_score) {
      max_taxa.clear();
      max_score = score;
      max_taxon = taxon;
      //LOG_DEBUG("max score: %u. max taxon: %u\n", max_score, max_taxon);
    } else if (score == max_score) {
      if (max_taxa.empty()) // Is this check needed?
        max_taxa.insert(max_taxon);
      max_taxa.insert(taxon);
    }
  }
  // If two LTR paths are tied for max, return LCA of all
  if(max_taxa.size()) {
    //LOG_DEBUG("Ambiguous. Get the lca of all of these. Size: %zu\n", max_taxa.size());
    auto sit(max_taxa.begin());
    max_taxon = *sit;
    for(++sit; sit != max_taxa.end(); ++sit) max_taxon = lca(parent_map, max_taxon, *sit);
  }

  return max_taxon;
}

std::string rand_string(std::size_t n) {
    std::string ret;
    ret.reserve(n);
    static const char set[] = "abcdefghijklmnopqrstuvwxyz123456";
    while(ret.size() < n) ret += set[std::rand() & 31];
    assert(ret.size() == n);
    return ret;
}

std::string get_firstline(const char *fn) {
    gzFile fp(gzopen(fn, "rb"));
    if(fp == nullptr) LOG_EXIT("Could not read from file %s\n", fn);
    static const std::size_t bufsz(2048);
    char buf[bufsz];
    std::string ret(gzgets(fp, buf, bufsz));
    if(ret.back() != '\n')
        throw std::runtime_error(std::string("[E:get_firstline] line from ") +
                                 fn + " did not fit in buffer of size " +
                                 std::to_string(bufsz) +
                                 ". Try recompiling with a larger "
                                 "buffer or rewrite get_firstline.");
    ret.pop_back();
    gzclose(fp);
    return ret;
}

tax_t get_taxid(const char *fn, khash_t(name) *name_hash) {
    gzFile fp(gzopen(fn, "rb"));
    if(fp == nullptr) LOG_EXIT("Could not read from file %s\n", fn);
    static const std::size_t bufsz(2048);
    khint_t ki;
    char buf[bufsz];
    char *line(gzgets(fp, buf, bufsz));
    char *p(++line);
    if(std::strchr(p, '|')) {
        p = std::strrchr(p, '|');
        while(*--p != '|');
        char *q(std::strchr(++p, '|'));
        *q = 0;
        line = p;
        //if(strchr(p, '.')) *strchr(p, '.') = 0;
    } else {
        while(!std::isspace(*p)) ++p;
        *p = 0;
    }
    if(unlikely((ki = kh_get(name, name_hash, line)) == kh_end(name_hash))) {
        LOG_WARNING("Missing taxid for %s, '%s'.\n", buf, line);
        gzclose(fp);
        return UINT32_C(-1);
    }
    LOG_DEBUG("Successfully got taxid %u for path %s\n", kh_val(name_hash, ki), fn);
    gzclose(fp);
    return kh_val(name_hash, ki);
}

std::map<uint32_t, uint32_t> kh2kr(khash_t(p) *map) {
    std::map<uint32_t, uint32_t> ret;
    if(map) 
        for(khiter_t ki(0); ki != kh_end(map); ++ki)
            if(kh_exist(map, ki))
                ret[kh_key(map, ki)] = kh_val(map, ki);
    return ret;
}



std::unordered_map<tax_t, strlist> tax2genome_map(khash_t(name) *name_map, const std::vector<std::string> &paths) {
    tax_t taxid;
    std::unordered_map<tax_t, std::forward_list<std::string>> ret;
    typename std::unordered_map<tax_t, std::forward_list<std::string>>::iterator m;
    ret.reserve(paths.size());
    for(const auto &path: paths) {
        taxid = get_taxid(path.data(), name_map);
        if(taxid == UINT32_C(-1)) continue;
        LOG_INFO("Found taxid %u\n", taxid);
        if((m = ret.find(taxid)) == ret.end()) ret.emplace(taxid, strlist{path});
        else                                   m->second.emplace_front(path);
    }
    return ret;
}

namespace {
template<typename T> class TD;
}


std::unordered_map<tax_t, std::set<tax_t>> make_ptc_map(const khash_t(p) *taxmap, const std::vector<tax_t> &taxes, const std::unordered_map<tax_t, ClassLevel> &lvl_map) {
    std::unordered_map<tax_t, std::set<tax_t>> ret;
    ret.reserve(kh_size(taxmap));
    std::vector<tax_t> sorted_taxes = taxes;
    SORT(sorted_taxes.begin(), sorted_taxes.end(), [&lvl_map](const tax_t a, const tax_t b) {
             return lvl_map.at(a) < lvl_map.at(b);
    });
#if !NDEBUG
    for(auto it1(sorted_taxes.begin()), it2(it1 + 1); it2 != sorted_taxes.end(); ++it1, ++it2) {
        assert(lvl_map.at(*it1) <= lvl_map.at(*it2));
    }
#endif
    typename std::unordered_map<tax_t, std::set<tax_t>>::iterator it;
    khiter_t ki;
    tax_t tmptax;
    for(const auto tax: sorted_taxes) {
#if !NDEBUG
        std::cerr << "Processing tax " << tax << '\n';
#endif
        if((it = ret.find(tax)) == ret.end()) ret.emplace(tax, std::set<tax_t>{});
        tmptax = tax;
        while((ki = kh_get(p, taxmap, tmptax)) != 0u) {
            if(ki == kh_end(taxmap)) {
                throw std::runtime_error(std::string("Missing parent for taxid ") + std::to_string(tax));
            }
            tmptax = kh_val(taxmap, ki);
#if !NDEBUG
            std::cerr << "Parent: " << tmptax << '\n';
#endif
            if((it = ret.find(tmptax)) == ret.end()) ret.emplace(tmptax, std::set<tax_t>{tax});
            else                                     it->second.insert(tax);
        }
    }
    return ret;
}

std::unordered_map<tax_t, strlist> tax2desc_genome_map(
        const std::unordered_map<tax_t, strlist> &tx2g,
        const khash_t(p) *taxmap, const std::vector<tax_t> &taxes,
        const std::unordered_map<tax_t, ClassLevel> &lvl_map) {
    std::unordered_map<tax_t, strlist> ret;
    for(const auto &pair: make_ptc_map(taxmap, taxes, lvl_map)) {
        typename std::unordered_map<tax_t, strlist>::const_iterator pit;
        strlist list;
        if((pit = tx2g.find(pair.first)) != tx2g.end()) list.insert_after(list.begin(), pit->second.begin(), pit->second.end());
        else {
            if(kh_get(p, taxmap, pair.first) == kh_end(taxmap)) {
                std::cerr << "No parent for node " << (int)pair.first << '\n';
                throw std::runtime_error(std::string("Invalid taxid ") + std::to_string(pair.first));
            } else {
#if !NDEBUG
                std::cerr << "Valid tax id " << kh_key(taxmap, kh_get(p, taxmap, pair.first))
                          << " with parent " << kh_val(taxmap, kh_get(p, taxmap, pair.first))
                          << ". Do nothing.\n";
#endif
            }
        }
        for(const auto child: pair.second)
            if((pit = tx2g.find(child)) != tx2g.end())
                list.insert_after(list.begin(), pit->second.begin(), pit->second.end());
        ret.emplace(pair.first, std::move(list));
    }
    return ret;
}

bool isfile(const char *path) noexcept {
    return access(path, F_OK) != -1;
}

const char *bool2str(bool val) {
    static const char * vals [] {
        "false",
        "true"
    };
    return vals[val];
}

#define _KHD(x) template<> void khash_destroy(khash_t(x) *map) noexcept {kh_destroy(x, map);}

_KHD(all)
_KHD(c)
_KHD(64)
_KHD(p)

ClassLevel get_linelvl(const char *line, std::string &buffer, const std::unordered_map<std::string, ClassLevel> &map) {
    const char *p(strchr(line, '|'));
    if(!p || (p = strchr(p + 1, '|')) == nullptr)
        throw std::runtime_error("Improperly formatted line");
    p = p + 2;
    const char *q(p);
    while(*q != '\t' && *q) ++q;
    if(!*q) throw std::runtime_error("Improperly formatted line");
    buffer = std::string(p, q);
    auto m(map.find(buffer));
    if(m == map.end()) {
        for(const auto &pair: map) {
            std::cerr << "Key: " << pair.first << ". Value: " << static_cast<int>(pair.second) << ".\n";
        }
        throw std::runtime_error(std::string("Unexpected field entry '") + buffer + "' for line " + line);
    }
    return m->second;
}

std::unordered_map<tax_t, ClassLevel> get_tax_depths(const khash_t(p) *taxmap, const char *path) {
    std::unordered_map<tax_t, ClassLevel> ret;
    std::ifstream ifs(path);
    std::string buffer;
    tax_t t;
    if(!ifs.good()) throw "a party";
    for(std::string line; std::getline(ifs, line);) {
        t = atoi(line.data());
        if(kh_get(p, taxmap, t) == kh_end(taxmap)) continue;
        ret.emplace(t, get_linelvl(line.data(), buffer, classlvl_map));
    }
    return ret;
}

std::vector<tax_t> get_sorted_taxes(const khash_t(p) *taxmap, const char *path) {
    std::vector<tax_t> taxes;
    {
        std::unordered_set<tax_t> taxset;
        for(khiter_t ki(0); ki != kh_end(taxmap); ++ki)
            if(kh_exist(taxmap, ki))
                taxset.insert(kh_key(taxmap, ki));
        taxes = std::vector<tax_t>(taxset.begin(), taxset.end());
    }
    std::unordered_map<tax_t, ClassLevel> taxclassmap(get_tax_depths(taxmap, path));
    typename std::unordered_map<tax_t, ClassLevel>::iterator ma, mb;
    SORT(taxes.begin(), taxes.end(), [&tcm=taxclassmap,tm=taxmap,&ma,&mb](const tax_t a, const tax_t b) {
        return (ma = tcm.find(a)) == tcm.end() ? (mb = tcm.find(b)) == tcm.end() ? a < b: false
                                               : (mb = tcm.find(b)) == tcm.end() ? true
                                                                                 : ma->second == mb->second ? get_parent(tm, a) < get_parent(tm, b)
                                                                                                            : ma->second > mb->second;
    });
    return taxes;
}

} //namespace emp
