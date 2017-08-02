#!/usr/bin/env python
import sys
import unittest
from array import array as pa  # Python array
import numpy as np


class TaxEntry(object):
    def __init__(self, id, parent, name=None, depth=-1)
        self.name = name
        self.id = id
        self.depth = depth
        self.parent = parent
        self.genome_paths = []

    def __hash__(self):
        return hash(self.id)


class Taxonomy(object):
    def __init__(self):
        self.taxes = set()


    def add_file(self, path):
        for line in open(path):
            add_line(line)

    def add_line(self, line):
        toks = [i.strip() for i in line.split('|')]
        self.taxes[int(toks[0])] = TaxEntry(int(toks[0]), int(toks[1]))

    def add_genome_path(self, path, tax):
        try:
            assert isinstance(tax, int)
            self.taxes[tax].genome_paths.append(path)
        except AssertionError:
            print("Assertion failed. Tax is not an int. Tax: %s, %s" % (str(tax), repr(tax)))
            raise
        except KeyError:
            print("Missing tax %i" % tax)
            raise


if __name__ == "__main__":
    class TC1(unittest.TestCase):

        def test_testing(self): pass
        def test_assert(self): self.assertTrue(True)
        def test_assert_false(self): self.assertFalse(False)

    unittest.main()
    raise NotImplementedError("Unit tests for this module are not complete.")
