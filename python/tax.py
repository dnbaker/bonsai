#!/usr/bin/env python
from six import iteritems
import sys
import unittest
from array import array as pa  # Python array
from . import kmer
import numpy as np
import itertools
from download_genomes import xfirstline


def load_nameidmap(path):
    with open(path) as f:
        return {el[0]: int(el[1]) for el in
                map(lambda x: x.split('\t'), f)}


class TaxEntry(object):
    def __init__(self, line, classlvl_map):
        toks = [i.strip() for i in line.split('|')]
        self.depth = classlvl_map[toks[2]]
        self.name = None
        self.id = int(toks[0])
        self.parent = int(toks[1])
        self.genome_paths = []

    def update(self, path):
        self.genome_paths.append(path)

    def __hash__(self):
        '''Hash function for TaxEntry.
           Note: this assumes that the user does not re-use tax ids.'''
        return hash(self.id)

    def __eq__(self, other):
        if isinstance(self, type(other)):
            return self.id == other.id
        raise NotImplementedError("Cannot compare TaxEntry "
                                  "objects with other objects.")

    def __cmp__(self, other):
        if isinstance(self, type(other)):
            if self.depth < other.depth:
                return 1
            if self.depth > other.depth:
                return -1
            return 0
        raise NotImplementedError("TaxEntry cannot be compared "
                                  "to %s" % str(type(other)))


class Taxonomy(object):

    def __init__(self, gi2taxpath, nodespath):
        self.taxes = []
        self.gi2taxmap = {line.split()[0]: int(line.split()[1]) for
                          line in open(gi2taxpath)}
        self.nodes_path = nodespath
        add_file(self.nodes_path)

    def add_file(self, path):
        for line in open(path):
            add_line(line)

    def get_parents(self, el):
        parent = self.__getitem__(el)
        ret = {parent}
        while parent:
            parent = self.__getitem__(el)
            ret.add(parent)
        return ret

    def build_child_ancestor_map():
        # This would take some careful thinking.
        raise NotImplementedError("Stuff")

    def __getitem__(self, el):
        return self.taxes[el]

    def __setitem__(self, key, val):
        self.taxes[key] = val

    def __contains__(self, key):
        return key in self.taxes

    def add_line(self, line):
        taxes.append(TaxEntry(line, self.gi2taxmap))

    def add_genome_path(self, path, tax):
        try:
            assert isinstance(tax, int)
            self.taxes[tax].genome_paths.append(path)
        except AssertionError:
            print("Assertion failed. Tax is not an int. Tax: %s, %s" %
                  (str(tax), repr(tax)))
            raise
        except KeyError:
            print("Missing tax %i" % tax)
            raise


if __name__ == "__main__":
    NAMEID_MAP_PATH = "save/combined_nip.txt"

    class TC1(unittest.TestCase):

        def test_testing(self):
            pass

        def test_assert(self):
            self.assertTrue(True)

        def test_assert_false(self):
            self.assertFalse(False)

        def test_load_nameid(self):
            map = load_nameidmap(NAMEID_MAP_PATH)
            for k, v in iteritems(map):
                self.assertIsInstance(v, int)
                self.assertIsInstance(k, str)

        def test_taxonomy_build(self):
            NODES_PATH = "ref/nodes.dmp"
            tax = Taxonomy()
            with open(NODES_PATH) as f:
                for line in f:
                    tax.add_line(line)
            self.assertTrue(0 in tax)

    unittest.main()
