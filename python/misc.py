#!/usr/bin/env python
import sys
import string
from collections import defaultdict


def freq(iterable):
    """
    Returns a dictionary of counts for each item in an iterable.
    >>>freq("ACGTTTAAA")
    {'A': 4, 'C': 1, 'G': 1, 'T': 3}
    """
    ret = defaultdict(int)
    for el in iterable:
        ret[el] += 1
    return ret

try:
    from cytoolz import frequencies as freq
except ImportError:
    pass
    # Don't sweat it

REV_CMP_TABLE = (str if sys.version_info[0] == 3
                 else string).maketrans("ACGTN", "TGCAN")


def revcmp(seq):
    """
    Returns the reverse complement of a sequence.
    >>>revcmp("ACGTNTTTAAATTT")
    'AAATTTAAANACGT'
    """
    return seq[::-1].translate(REV_CMP_TABLE)


def xopen(path):
    """
        Stolen from Dooplicity. (https://github.com/nellore/rail/),
        then stripped to only open files with open or gzip to open
        based on magic number presence.
    """
    import gzip
    fh = (gzip.open(path, "rb") if open(path, 'rb').read(2) == '\x1f\x8b'
          else open(path, "r"))
    try:
        yield fh
    finally:
        fh.close()

__all__ = [revcmp, REV_CMP_TABLE, freq, xopen]


if __name__ == "__main__":
    """
    Unit tests
    """
    import unittest

    class Test(unittest.TestCase):

        def test_revcmp(self):
            self.assertEqual(revcmp("ACGTACCTTATATATATA"),
                             "TATATATATAAGGTACGT")

        def test_freq(self):
            self.assertEqual(freq("ACGTTTAAA"),
                             {'A': 4, 'C': 1, 'G': 1, 'T': 3})
    unittest.main()
