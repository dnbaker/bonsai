from __future__ import print_function
from array import array
import pysam
import sys

KMER_LUT = array('B', [255] * 256)
KMER_LUT[ord('A')] = 0
KMER_LUT[ord('C')] = 1
KMER_LUT[ord('G')] = 2
KMER_LUT[ord('T')] = 3
KMER_LUT[ord('a')] = 0
KMER_LUT[ord('c')] = 1
KMER_LUT[ord('g')] = 2
KMER_LUT[ord('t')] = 3
AMBIGUOUS = -1  # Python has the sign bit. Use it.


def str2kmerint(s):
    ret = 0
    for i in range(len(s)):
        ret <<= 2
        val = KMER_LUT[ord(s[i])]
        if val == 255:
            return AMBIGUOUS
        ret |= val
    return ret


def kmer2str(kmer, k=-1):
    ret = ""
    while kmer:
        ret += "ACGT"[kmer & 3]
        kmer >>= 2
    return ret[::-1] if k < 0 else 'A' * (k - len(ret)) + ret[::-1]


def genome2kmergen(path, k=31):
    def __gen_seqs(name, k):
        for record in pysam.FastxFile(name):
            seq = record.sequence
            for i in range(len(seq) - k + 1):
                yield str2kmerint(seq[i:i + k])
    return __gen_seqs(path, k)


def genome2kmerset(path, k=31):
    # Will this get more expensive on Python 2 bc map?
    ret = set(genome2kmergen(path, k))
    if AMBIGUOUS in ret:
        ret.remove(AMBIGUOUS)
    return ret


__all__ = [KMER_LUT, AMBIGUOUS, str2kmerint, genome2kmergen,
           genome2kmerset, kmer2str]


if __name__ == "__main__":
    import unittest

    class UnitTests(unittest.TestCase):
        def setUp(self):
            pass

        def test_stuff(self):
            self.assertFalse(False)

        def test_kmergen(self):
            k31 = genome2kmerset("test/phix.fa", 31)
            print(k31)
            k180 = genome2kmerset("test/phix.fa", 180)
            print(k180)
            k1024 = genome2kmerset("test/phix.fa", 1024)
            print(k1024)
            for name, item in zip(("k31", "k180", "k1024"),
                                  (k31, k180, k1024)):
                print("Length of %s is %i\n" % (name, len(item)))

    unittest.main()
