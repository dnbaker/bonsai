#include "util.h"
#include <set>
#include <cassert>

namespace emp {

std::size_t count_lines(const char *fn) noexcept {
    std::FILE *fp(std::fopen(fn, "r"));
    std::size_t bufsz = 4096;
    char *buf((char *)std::malloc(bufsz));
    ssize_t len;
    std::size_t n(0);
    while((len = getline(&buf, &bufsz, fp)) >= 0) ++n;
    free(buf);
    std::fclose(fp);
    return n;
}

// Rewritten from Kraken's source code.
// Consider rewriting to use kbtree instead of std::map.
// lh3's benchmarks indicate that his is only as fast as std::map,
// though it uses less than half the memory. The thing is, taxonomic trees
// aren't that deep.
// We don't gain much from optimizing that.
std::uint32_t lca(khash_t(p) *map, std::uint32_t a, std::uint32_t b) noexcept {
    // Use std::set to use RB trees for small set rather than hash table.
    std::set<std::uint32_t> nodes;
    if(a == b) return a;
    khint_t ki;
    while(a) {
        nodes.insert(a);
        if((ki = kh_get(p, map, a)) == kh_end(map)) {
            return (std::uint32_t)-1;
        }
        a = kh_val(map, ki);
    }
    while(b) {
        if(nodes.find(b) != nodes.end()) {
            return b;
        }
        if((ki = kh_get(p, map, b)) == kh_end(map)) {
            fprintf(stderr, "Missing taxid %u. Returning -1!\n", b);
            return (std::uint32_t)-1;
        }
        b = kh_val(map, ki);
    }
    return 1;
}

unsigned node_depth(khash_t(p) *map, std::uint32_t a) noexcept {
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
    free(buf);
    return ret;
}

void destroy_name_hash(khash_t(name) *hash) noexcept {
    for(khint_t ki(kh_begin(hash)); ki != kh_end(hash); ++ki)
        if(kh_exist(hash, ki))
            free((void *)kh_key(hash, ki));
    kh_destroy(name, hash);
}

khash_t(p) *build_parent_map(const char *fn) noexcept {
    std::size_t nlines(count_lines(fn));
    std::FILE *fp(std::fopen(fn, "r"));
    khash_t(p) *ret(kh_init(p));
    kh_resize(p, ret, nlines);
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
    free(buf);
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
std::uint32_t resolve_tree(std::map<std::uint32_t, std::uint32_t> &hit_counts,
                      khash_t(p) *parent_map) noexcept
{
  std::set<std::uint32_t> max_taxa;
  std::uint32_t max_taxon(0), max_score(0);
  auto it(hit_counts.begin());

  // Sum each taxon's LTR path
  while (it != hit_counts.end()) {
    std::uint32_t taxon(it->first), node(taxon), score(0);
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
    ++it;
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


#define _KHD(x) template<> void khash_destroy(khash_t(x) *map) noexcept {kh_destroy(x, map);}

_KHD(all)
_KHD(c)
_KHD(64)
_KHD(p)

} //namespace emp
