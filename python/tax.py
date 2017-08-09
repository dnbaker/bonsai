#!/usr/bin/env python
from six import iteritems
import sys
import unittest
from array import array as pa  # Python array
import kmer
import numpy as np
import itertools
from download_genomes import xfirstline
from parse_nodes import generate_python_class_map


def load_nameidmap(path):
    with open(path) as f:
        return {el[0]: int(el[1]) for el in
                map(lambda x: x.split('\t'), f)}


class TaxEntry(object):
    def __init__(self, line, gi2taxmap, classlvl_map,
                 __val_only1=None, __val_only2=None):
        if isinstance(line, int):
            self.id = line
            self.parent = gi2taxmap
            self.depth = classlvl_map
            self.genome_paths = (__val_only1 if isinstance(__val_only1, list)
                                 else None)
            self.name = __val_only2
            if self.id == 1:
                self.parent = 0
        elif isinstance(line, str) is False:
            print("Type of line is '%s'" % type(line))
            raise RuntimeError("ZOMGZ")
        else:
            toks = (i.strip() for i in line.split('|'))
            self.id = int(next(toks))
            self.parent = int(next(toks))
            try:
                self.depth = classlvl_map[next(toks)]
            except KeyError:
                print("'" + "', '".join(str(i) for
                      i in classlvl_map.keys()) + "'")
                raise
            self.name = None
            self.genome_paths = []

    def __str__(self):
        return "TaxEntry(%i, %i, %i, [%s], %s)" % (
            self.id, self.parent, self.depth,
            ", ".join("'%s'" % path for path in self.genome_paths),
            '"%s"' % self.name if self.name else "None")

    def update(self, path):
        self.genome_paths.append(path)

    def __hash__(self):
        '''Hash function for TaxEntry.
           Note: this assumes that the user does not re-use tax ids.'''
        return hash(self.id)

    def __eq__(self, other):
        return (isinstance(self, type(other)) and self.id == other.id and
                self.parent == other.parent)

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
        self.index = 0
        self.lvl_map = generate_python_class_map()
        try:
            assert "superkingdom" in self.lvl_map
        except AssertionError:
            print("Map string keys: " + ", ".join(i for i in
                                                  self.lvl_map.keys() if
                                                  isinstance(i, str)))
            print("Map int keys: " + ", ".join(str(i) for i in
                                               self.lvl_map.keys() if
                                               isinstance(i, int)))
            raise
        self.taxes = [TaxEntry(0, 0, 0, [], None)]
        with open(gi2taxpath) as f:
            self.gi2taxmap = {line.split()[0]: int(line.split()[1]) for
                              line in f}
        self.nodes_path = nodespath
        self.add_file(self.nodes_path)

    def __iter__(self):
        return self

    def __next__(self):
        try:
            ret = self.taxes[self.index]
            self.index += 1
        except IndexError:
            self.index = 0
            raise StopIteration()
        return ret

    def add_file(self, path):
        with open(path) as f:
            for line in f:
                self.add_line(line)

    def get_parents(self, el):
        parent = self.__getitem__(el)
        ret = {parent}
        while parent:
            parent = self.__getitem__(el)
            ret.add(parent)
        return ret

    def build_child_ancestor_map():
        ret = {key: [] for key in self.taxes}
        self.taxes.sort()
        print("Sorted taxes! Now:\n\n\n")
        print(", ".join(map(str, self.taxes)) + "\n\n\n")
        # This will take some careful thinking.
        raise NotImplementedError("Stuff. But save sorted taxes in this "
                                  "exception trace: %s" %
                                  ", ".join(map(str, self.taxes)) + "\n\n\n")
        return ret

    def __getitem__(self, el):
        try:
            return next(x for x in self.taxes if x == el or x.id == el)
        except StopIteration:
            raise KeyError("Missing element or id %s" % str(el))

    def __setitem__(self, key, val):
        try:
            next(x for x in self.taxes if x == el or
                 x.id == el).genome_paths.append(val)
        except StopIteration:
            raise KeyError("Missing element or id %s" % str(el))

    def __contains__(self, key):
        return (key in self.taxes if isinstance(key, TaxEntry) else
                any(tax.id == key for tax in self.taxes))

    def add_line(self, line):
        self.taxes.append(TaxEntry(line, self.gi2taxmap, self.lvl_map))

    def add_genome_path(self, path, tax):
        try:
            assert isinstance(tax, int)
            next(x for x in self.taxes if
                 x.id == tax).genome_paths.append(path)
        except AssertionError:
            print("Assertion failed. Tax is not an int. Tax: %s, %s" %
                  (str(tax), repr(tax)))
            raise
        except StopIteration:
            print("Missing tax %i" % tax)
            raise


if __name__ == "__main__":
    NAMEID_MAP_PATH = "save/combined_nip.txt"
    NODES_PATH = "ref/nodes.dmp"

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
            tax = Taxonomy(NAMEID_MAP_PATH, NODES_PATH)
            for el in tax:
                print("Element %s in tax" % el)
                el2 = eval("%s" % el)
                assert el == el2
            els = set(tax)
            print("els: '%s'" % (", ".join(map(str, els))))
            if 0 not in tax:
                print("Failed!")
                print(", ".join(i for i in tax))
            self.assertTrue(0 in tax)
            tax.build_child_ancestor_map()

    unittest.main()
