#include "util.h"
#include <set>
#include <cassert>
#include <cstring>
#include <zlib.h>
#include <sstream>
#include <fstream>
#include "kspp/ks.h"

namespace emp {

size_t count_lines(const char *fn) noexcept {
    std::ifstream is(fn);
    if(!is.good()) LOG_EXIT("Could not open file at %s\n", fn);
    std::string line;
    size_t n(0);
    while(std::getline(is, line)) ++n;
    return n;
}

std::unordered_map<tax_t, std::vector<tax_t>> invert_parent_map(const khash_t(p) *tax) noexcept {
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
    if(unlikely(map == nullptr)) {
        std::fprintf(stderr, "null taxonomy.\n");
        std::exit(EXIT_FAILURE);
    }
    std::set<tax_t> nodes;
    if(a == b) return a;
    if(b == 0) return a;
    if(a == 0) return b;
    khint_t ki;
    while(a) {
        nodes.insert(a);
        if((ki = kh_get(p, map, a)) == kh_end(map)) {
            std::fprintf(stderr, "Missing taxid %u. Returning node b (%u)!\n", a, b);
            return (tax_t)-1;
        }
        a = kh_val(map, ki);
    }
    while(b) {
        if(nodes.find(b) != nodes.end()) return b;
        if((ki = kh_get(p, map, b)) == kh_end(map)) {
            std::fprintf(stderr, "Missing taxid %u. Returning node b (%u)!\n", b, a);
            return (tax_t)-1;
        }
        b = kh_val(map, ki);
    }
    return 1;
}

unsigned node_dist(const khash_t(p) *map, tax_t leaf, tax_t root) noexcept {
    unsigned ret(0);
    khint_t ki;
    while(leaf) {
        if((ki = kh_get(p, map, leaf)) == kh_end(map)) {
            std::fprintf(stderr, "Tax ID %u missing. Abort!\n", leaf);
            std::exit(1);
        }
        ++ret;
        if((leaf = kh_val(map, ki)) == root) return ret;
    }
    LOG_EXIT("leaf %u is not a child of root %u\n", leaf, root);
    return ret;
}
unsigned node_depth(const khash_t(p) *map, tax_t a) noexcept {
    unsigned ret(0);
    khint_t ki;
    while(a) {
        if((ki = kh_get(p, map, a)) == kh_end(map)) {
            std::fprintf(stderr, "Tax ID %u missing. Abort!\n", a);
            std::exit(1);
        }
        a = kh_val(map, ki);
        ++ret;
    }
    return ret;
}

khash_t(name) *build_name_hash(const char *fn) noexcept {
    size_t bufsz(2048), namelen;
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
        p = ::emp::strchrnul(buf, '\t');
        ki = kh_put(name, ret, buf, &khr);
        if(khr == 0) { // Key already present.
            LOG_INFO("Key %s already present. Updating value from "
                     "%i to %s.\tNote: if you have performed TaxonomyReformation, this is an error.\n", kh_key(ret, ki), kh_val(ret, ki), p + 1);
        } else {
            namelen = p - buf;
            char *tmp = static_cast<char *>(std::malloc(namelen + 1));
            std::memcpy(tmp, buf, namelen);
            tmp[namelen] = '\0';
            kh_key(ret, ki) = tmp;
        }
        kh_val(ret, ki) = std::atoi(p + 1);
    }
    std::fclose(fp);
    std::free(buf);
    return ret;
}

void destroy_name_hash(khash_t(name) *hash) noexcept {
    if(hash == nullptr) return;
    for(khint_t ki(kh_begin(hash)); ki != kh_end(hash); ++ki)
        if(kh_exist(hash, ki))
            std::free((void *)kh_key(hash, ki));
    kh_destroy(name, hash);
}

std::map<tax_t, tax_t> build_kraken_tax(const std::string &fname) {
    const char *fn(fname.data());
    std::FILE *fp(std::fopen(fn, "r"));
    std::map<tax_t, tax_t> ret;
    size_t bufsz = 4096;
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

tax_t lca(const std::map<tax_t, tax_t> &parent_map, tax_t a, tax_t b)
{
  std::map<tax_t, tax_t>::const_iterator m;
  if(a == 0 || b == 0) return a ? a : b;

  std::set<uint32_t> a_path;
  while (a) a_path.insert(a), a = parent_map.find(a)->second;
  while (b) {
    if(a_path.find(b) != a_path.end())
      return b;
    b = parent_map.find(b)->second;
  }
  return 1;
}

khash_t(p) *build_parent_map(const char *fn) noexcept {
    std::ifstream is(fn);
    khash_t(p) *ret(kh_init(p));
    khint_t ki;
    int khr;
    std::string line;
    char *p;
    while(std::getline(is, line)) {
        switch(line[0]) case '\n': case '0': case '#': continue;
        ki = kh_put(p, ret, std::atoi(line.data()), &khr);
        kh_val(ret, ki) = (p = std::strchr(line.data(), '|')) ? std::atoi(p + 2)
                                                              : tax_t(-1);
        if(kh_val(ret, ki) == tax_t(-1)) LOG_WARNING("Malformed line: %s", line.data());
    }
    ki = kh_put(p, ret, 1, &khr);
    kh_val(ret, ki) = 0; // Root of the tree.
    LOG_DEBUG("Built parent map of size %zu from path %s\n", kh_size(ret), fn);
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
 // Taken, modified from Kraken with new, faster containers.
tax_t resolve_tree(const std::map<tax_t, tax_t> &hit_counts,
                   const khash_t(p) *parent_map) noexcept
{
  linear::set<tax_t> max_taxa;
  tax_t max_taxon(0), max_score(0);

  // Sum each taxon's LTR path
  for(auto it(hit_counts.cbegin()), e(hit_counts.cend()); it != e; ++it) {
    tax_t taxon(it->first), node(taxon), score(0);
    // Instead of while node > 0
    while(node)
        score += hit_counts.at(node), node = kh_val(parent_map, kh_get(p, parent_map, node));
    if(score > max_score) {
      max_taxa.clear();
      max_score = score;
      max_taxon = taxon;
      //LOG_DEBUG("max score: %u. max taxon: %u\n", max_score, max_taxon);
    } else if(score == max_score) {
      if(max_taxa.empty()) // Is this check needed?
        max_taxa.insert(max_taxon);
      max_taxa.insert(taxon);
    }
  }
  // If two LTR paths are tied for max, return LCA of all
  if(max_taxa.size()) {
    auto sit(max_taxa.begin());
    for(max_taxon = *sit++;sit != max_taxa.end(); max_taxon = lca(parent_map, max_taxon, *sit++));
  }

  return max_taxon;
}

tax_t resolve_tree(const linear::counter<tax_t, u16> &hit_counts,
                   const khash_t(p) *parent_map) noexcept
{
  linear::set<tax_t> max_taxa;
  tax_t max_taxon(0), max_score(0);

  // Sum each taxon's LTR path
  for(unsigned i(0); i < hit_counts.size(); ++i) {
    tax_t taxon(hit_counts.keys()[i]), node(taxon), score(0);
    // Instead of while node > 0
    while(node) {
        score += hit_counts.count(node);
        node = kh_val(parent_map, kh_get(p, parent_map, node));
    }
    // In C++20, we can use the spaceship operator!
    if(score > max_score) {
      max_taxa.clear();
      max_score = score;
      max_taxon = taxon;
    } else if(score == max_score) {
      if(max_taxa.empty()) // Is this check needed?
        max_taxa.insert(max_taxon);
      max_taxa.insert(taxon);
    }
  }
  // If two LTR paths are tied for max, return LCA of all
  if(max_taxa.size()) {
    auto sit(max_taxa.begin());
    for(max_taxon = *sit++;sit != max_taxa.end(); max_taxon = lca(parent_map, max_taxon, *sit++));
  }

  return max_taxon;
}


std::string rand_string(size_t n) {
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
    static const size_t bufsz(2048);
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

tax_t get_taxid(const char *fn, const khash_t(name) *name_hash) {
    gzFile fp(gzopen(fn, "rb"));
    if(fp == nullptr) LOG_EXIT("Could not read from file %s\n", fn);
    static const size_t bufsz(2048);
    khint_t ki;
    char buf[bufsz];
    char *line(gzgets(fp, buf, bufsz));
    if(line == nullptr) {
        int err;
        LOG_INFO("zlib error: %s\n", gzerror(fp, &err));
        throw zlib_error(err, fn);
    }
    char *p(++line);
#ifdef SYNTHETIC_GENOME_EXPERIMENTS
    const tax_t ret(std::atoi(p));
#else
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
    const tax_t ret(unlikely((ki = kh_get(name, name_hash, line)) == kh_end(name_hash)) ? UINT32_C(-1) : kh_val(name_hash, ki));
#endif
    gzclose(fp);
    return ret;
}

std::map<uint32_t, uint32_t> kh2kr(khash_t(p) *map) {
    std::map<uint32_t, uint32_t> ret;
    if(map)
        for(khiter_t ki(0); ki != kh_end(map); ++ki)
            if(kh_exist(map, ki))
                ret[kh_key(map, ki)] = kh_val(map, ki);
    return ret;
}
#if 0
#define PRINT_LISTMAP(ret) do \
    {\
        for(const auto &pair: ret) {\
            std::cerr << "key: " << pair.first << ',' << "children: ";\
            for(const auto el: pair.second)\
                std::string tmp = el;\
                std::cerr << pair.first << ',';\
            std::cerr << '\n';\
        }\
    } while(0)

#else
#define PRINT_LISTMAP(ret)
#endif


std::unordered_map<tax_t, strlist> tax2genome_map(khash_t(name) *name_map, const std::vector<std::string> &paths) {
    tax_t taxid;
    std::unordered_map<tax_t, std::forward_list<std::string>> ret;
    typename std::unordered_map<tax_t, std::forward_list<std::string>>::iterator m;
    ret.reserve(paths.size());
#if !NDEBUG
    ks::string ks;
#endif
    for(const auto &path: paths) {
        if((taxid = get_taxid(path.data(), name_map)) == UINT32_C(-1)) continue;
        if((m = ret.find(taxid)) == ret.end()) m = ret.emplace(taxid, strlist{path}).first;
        else if(std::find(m->second.begin(), m->second.end(), path) == m->second.end()) m->second.push_front(path);
#if !NDEBUG
        ks.clear();
        ks.sprintf("Path %s has taxid %u.\n", path.data(), taxid);
        for(const std::string &p: m->second) {
            ks.sprintf("Element: \"%s\",", p.data());
        }
        ks.back() = '\n';
        std::cerr << ks.data();
        if((ret.size() & (ret.size() - 1)) == 0) std::cerr << "Size of tax2genome_map is now " << ret.size() << '\n';
        PRINT_LISTMAP(ret);
#endif
    }
    PRINT_LISTMAP(ret);
    return ret;
}


std::unordered_map<tax_t, std::set<tax_t>> make_ptc_map(
        const khash_t(p) *taxmap, const std::vector<tax_t> &taxes,
        const std::unordered_map<tax_t, ClassLevel> &lvl_map) {
    std::unordered_map<tax_t, std::set<tax_t>> ret;
    ret.reserve(kh_size(taxmap));
    std::vector<tax_t> sorted_taxes = taxes;
#if !NDEBUG
    bool fail(false);
    for(const auto tax: sorted_taxes) {
        if(lvl_map.find(tax) == lvl_map.end()) std::cerr << "Missing tax level for " << tax << '\n', fail = true;
    }
    if(fail) throw std::runtime_error("Failed for missin tax levels.");
#endif
    ::std::sort(sorted_taxes.begin(), sorted_taxes.end(), [&lvl_map](const tax_t a, const tax_t b) {
             try {
                return lvl_map.at(a) < lvl_map.at(b);
             } catch(std::out_of_range &ex) {
                std::cerr << "Out of range: " << ex.what() << '\n';
                if(lvl_map.find(a) == lvl_map.end()) {std::cerr << "Missing tax " << a << '\n'; throw;}
                if(lvl_map.find(b) == lvl_map.end()) {std::cerr << "Missing tax " << b << '\n'; throw;}
                throw std::runtime_error("ZOMG");
             }
    });
#if !NDEBUG
    for(auto it1(sorted_taxes.begin()), it2(it1 + 1); it2 != sorted_taxes.end(); ++it1, ++it2) {
        assert(lvl_map.at(*it1) <= lvl_map.at(*it2));
    }
#endif
    typename std::unordered_map<tax_t, std::set<tax_t>>::iterator it;
    khiter_t ki;
    tax_t tmptax;
    for(auto sit(sorted_taxes.begin()), eit(sorted_taxes.end()); sit != eit; ++sit) {
        const tax_t tax(*sit);
#if 0
        std::cerr << "Processing tax " << tax << '#' << static_cast<size_t>(sit - sorted_taxes.begin()) << " of " << sorted_taxes.size() << '\n';
#endif
        if((it = ret.find(tax)) == ret.end()) ret.emplace(tax, std::set<tax_t>{});
        tmptax = tax;
        while((ki = kh_get(p, taxmap, tmptax)) != kh_end(taxmap)) {
            if(kh_val(taxmap, ki) == 0) goto loop_end;
            tmptax = kh_val(taxmap, ki);
            if((it = ret.find(tmptax)) == ret.end()) ret.emplace(tmptax, std::set<tax_t>{tax});
            else                                     it->second.insert(tax);
        }
        // Only reach if taxid is missing from taxonomy file.
        throw std::runtime_error(std::string("Missing parent for taxid ") + std::to_string(tax));
        loop_end:
#if 0
        std::cerr << "Now finishing up for tax " << tax << '\n';
#else
        continue; // Syntactic boilerplate.
#endif
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
        if((pit = tx2g.find(pair.first)) != tx2g.end()) {
            for(const auto &str: pit->second) list.push_front(str);
        }
        else {
            if(kh_get(p, taxmap, pair.first) == kh_end(taxmap)) {
                std::cerr << "No parent for node " << (int)pair.first << '\n';
                throw std::runtime_error(std::string("Invalid taxid ") + std::to_string(pair.first));
            } else {
#if 0
                std::cerr << "Valid tax id " << kh_key(taxmap, kh_get(p, taxmap, pair.first))
                          << " with parent " << kh_val(taxmap, kh_get(p, taxmap, pair.first))
                          << ". Do nothing.\n";
#endif
            }
        }
        for(const auto child: pair.second)
            if((pit = tx2g.find(child)) != tx2g.end())
                for(const auto &el: pit->second)
                    list.push_front(el);
        ret.emplace(pair.first, std::move(list));
    }
    return ret;
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

template<> void khash_destroy(khash_t(name) *map) noexcept {destroy_name_hash(map);}

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
    if(!ifs.good()) throw std::runtime_error(std::string("could not open file at ") + path);
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

std::vector<tax_t> get_desc_lca(tax_t a, tax_t b, const std::unordered_map<tax_t, std::vector<tax_t>> &parent_map, const khash_t(p) *taxmap) {
    return get_all_descendents(parent_map, lca(taxmap, a, b));
}
void print_name_hash(khash_t(name) *hash) noexcept {
    for(khiter_t ki(0); ki < kh_size(hash); ++ki) {
        if(kh_exist(hash, ki)) {
            std::fprintf(stderr, "Key: %s. value: %u\n", kh_key(hash, ki), kh_val(hash, ki));
        }
    }
}

tax_t get_max_val(const khash_t(p) *hash) noexcept {
    tax_t mx(std::numeric_limits<tax_t>::min());
    for(khiter_t ki(0); ki < kh_size(hash); ++ki)
        if(kh_exist(hash, ki))
            mx = std::max(std::max(kh_key(hash, ki), kh_val(hash, ki)), mx);
    return mx;
}

khash_t(all) *load_binary_kmerset(const char *path) {
    std::FILE *fp(std::fopen(path, "rb"));
    if(fp == nullptr) throw std::system_error(std::error_code(2, std::system_category()), std::string("Cannot open path at ") + path + ".\n");
    khash_t(all) *ret(kh_init(all));
    u64 n;
    std::fread(&n, 1, sizeof(n), fp);
    if(kh_resize(all, ret, n) < 0) LOG_EXIT("Could not resize hash table to next power of 2 above %zu. New size: %zu\n", n, kh_n_buckets(ret));
    LOG_DEBUG("About to place %zu elements into a hash table of max size %zu\n", n, kh_n_buckets(ret));
    for(int khr; std::fread(&n, 1, sizeof(u64), fp) == sizeof(u64); kh_put(all, ret, n, &khr));
    std::fclose(fp);
#if !NDEBUG
    // Just make sure it all worked.
    for(khiter_t ki(0); ki < kh_end(ret); ++ki) {
        if(kh_exist(ret, ki)) assert(kh_get(all, ret, kh_key(ret, ki)) != kh_end(ret));
    }
    fp = std::fopen(path, "rb");
    std::fread(&n, 1, sizeof(u64), fp); // Skip first number.
    while(std::fread(&n, 1, sizeof(u64), fp) == sizeof(u64)) assert(kh_get(all, ret, n) != kh_end(ret));
    std::fclose(fp);
#endif
    return ret;
}


lazy::vector<u64, size_t> load_binary_kmers(const char *path) {
    std::FILE *fp(std::fopen(path, "rb"));
    if(fp == nullptr) throw std::system_error(std::error_code(2, std::system_category()), std::string("Cannot open path at ") + path + ".\n");
    lazy::vector<u64, size_t> ret;
    u64 n;
    std::fread(&n, sizeof(n), 1, fp);
    ret.resize(n, lazy::LAZY_VEC_NOINIT);
    auto it(ret.begin());
    auto eit(ret.end());
    while(std::fread(it++, sizeof(u64), 1, fp) == sizeof(u64))
        if(unlikely(it == eit))
            throw std::runtime_error("Read too many integers from file. Number provided (" + std::to_string(n) + ") is wrong.");
    std::fclose(fp);
    return ret;
}

std::vector<std::string> get_paths(const char *path) {
    std::ifstream is(path);
    std::vector<std::string> ret;
    for(std::string line;std::getline(is, line);ret.emplace_back(std::move(line)));
    return ret;
}

std::ifstream::pos_type filesize(const char* filename)
{
    std::ifstream in(filename, std::ifstream::ate | std::ifstream::binary);
    return in.tellg();
}

} //namespace emp
