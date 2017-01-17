#ifndef CUCKOO_FILTER_HASHUTIL_H_
#define CUCKOO_FILTER_HASHUTIL_H_

/* From https://github.com/efficient/cuckoofilter under the Apache license.
 *
*/

#include <stdint.h>
#include <stdlib.h>
#include <sys/types.h>
#include <cinttypes>

#include <string>

#include <openssl/evp.h>

namespace cuckoofilter {

/*      UINT64_MAX 18446744073709551615ULL */
#define P10_UINT64 10000000000000000000ULL /* 19 zeroes */
#define E10_UINT64 19

#define STRINGIZER(x) # x
#define TO_STRING(x) STRINGIZER(x)

static inline int putu128(__uint128_t big, FILE *fp) {
  size_t rc = 0;
  size_t i = 0;
  if (big >> 64) {
    char buf[40];
    while (big / P10_UINT64) {
      rc += sprintf(buf + E10_UINT64 * i, "%." TO_STRING(E10_UINT64) PRIu64, (uint64_t)(big % P10_UINT64));
      ++i;
      big /= P10_UINT64;
    }
    rc += fprintf(fp, "%" PRIu64, (uint64_t)big);
    while (i--) {
      fwrite(buf + E10_UINT64 * i, sizeof(char), E10_UINT64, fp);
      ++rc;
    }
  } else rc += fprintf(fp, "%" PRIu64, (uint64_t)big);
  return rc;
}

class HashUtil {
 public:
  // Bob Jenkins Hash
  static uint32_t BobHash(const void *buf, size_t length, uint32_t seed = 0);
  static uint32_t BobHash(const std::string &s, uint32_t seed = 0);

  // Bob Jenkins Hash that returns two indices in one call
  // Useful for Cuckoo hashing, power of two choices, etc.
  // Use idx1 before idx2, when possible. idx1 and idx2 should be initialized to
  // seeds.
  static void BobHash(const void *buf, size_t length, uint32_t *idx1,
                      uint32_t *idx2);
  static void BobHash(const std::string &s, uint32_t *idx1, uint32_t *idx2);

  // MurmurHash2
  static uint32_t MurmurHash(const void *buf, size_t length, uint32_t seed = 0);
  static uint32_t MurmurHash(const std::string &s, uint32_t seed = 0);

  // SuperFastHash
  static uint32_t SuperFastHash(const void *buf, size_t len);
  static uint32_t SuperFastHash(const std::string &s);

  // Null hash (shift and mask)
  static uint32_t NullHash(const void *buf, size_t length, uint32_t shiftbytes);

  // Wrappers for MD5 and SHA1 hashing using EVP
  static std::string MD5Hash(const char *inbuf, size_t in_length);
  static std::string SHA1Hash(const char *inbuf, size_t in_length);


  // See Martin Dietzfelbinger, "Universal hashing and k-wise independent random
  // variables via integer arithmetic without primes".
  static uint64_t TwoIndependentMultiplyShift(uint64_t key) {
    const uint64_t SEED[4] = {0x818c3f78ull, 0x672f4a3aull, 0xabd04d69ull,
                              0x12b51f95ull};
    const unsigned __int128 m =
        *reinterpret_cast<const unsigned __int128 *>(&SEED[0]);
/*
    const unsigned __int128 m =
        *reinterpret_cast<const unsigned __int128 *>(&SEED[0]);
*/
    const unsigned __int128 a =
        *reinterpret_cast<const unsigned __int128 *>(&SEED[2]);
    return (a + m * key) >> 64;
  }

 private:
  HashUtil();
};

}  // namespace cuckoofilter

#endif  // CUCKOO_FILTER_HASHUTIL_H_
