#!/usr/bin/env python
import sys
import unittest
import numpy as np
import itertools
import operator
from six import iteritems

import parse_nodes
import kmer

def load_nameidmap(path):
    with open(path) as f:
        return {el[0]: int(el[1]) for el in
                map(lambda x: x.split('\t'), f)}


class TaxEntry(object):
    def __init__(self, line, gi2taxmap, classlvl_map,
                 __val_only1=None, __val_only2=None):

        if isinstance(line, int):
            #  Only for creating the ROOT node. This is horrible.
            assert isinstance(classlvl_map, int)  # This should provide the taxonomic depth
            self.id = line
            self.parent = gi2taxmap
            self.depth = classlvl_map
            self.genome_paths = (__val_only1 if isinstance(__val_only1, list)
                                 else None)
            self.name = __val_only2
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
                print("key ('%s') missing from dictionary. Abort!", file=sys.stderr)
                raise
            self.name = None
            self.genome_paths = []
        if self.id == 1:
            self.parent = 0

    def __str__(self):
        return "TaxEntry(%i, %i, %i, [%s], %s)" % (
            self.id, self.parent, self.depth,
            ", ".join("'%s'" % path for path in self.genome_paths),
            '"%s"' % self.name if self.name else "None")

    def update(self, path):
        self.genome_paths.append(path)

    def __hash__(self):
        '''
        Hash function for TaxEntry.
        Multiple genomes for one taxid should be in the list "genome_paths".
        There should not be multiple TaxEntry objects for one taxid
        '''
        return hash(self.id)

    def __eq__(self, other):
        return (isinstance(self, type(other)) and self.id == other.id and
                self.parent == other.parent)

    def __lt__(self, other):
        if isinstance(self, type(other)):
            if self.depth > other.depth:
                return 1
            elif self.depth < other.depth:
                return 0
            if self.parent < other.parent:
                return 1
            elif self.parent > other.parent:
                return 0
            else:
                return self.id < other.id
        else:
            raise NotImplementedError("Cannot currently compare TaxEntry "
                                      "to %s" % type(other))

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
        self.lvl_map = parse_nodes.generate_python_class_map()
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
        sys.stderr.write("I successfully asserted that lvl_map belongs in existence.\n")
        self.taxes = [TaxEntry(0, 0, 0, [], "root")]
        with open(gi2taxpath) as f:
            self.gi2taxmap = dict(tuple(line.split()[:2]) for line in f)
        self.nodes_path = nodespath
        self.add_file(self.nodes_path)

    def __len__(self):
        return len(self.taxes)

    def __iter__(self):
        for i in self.taxes:
            yield i

    def __next__(self):
        for i in self.taxes:
            yield i
    '''
        try:
            ret = self.taxes[self.index]
            self.index += 1
        except IndexError:
            self.index = 0
            raise StopIteration()
        return ret
    '''

    def add_file(self, path):
        with open(path) as f:
            for line in f:
                self.add_line(line)

    def get_parents(self, el):
        parent = self.__getitem__(el)
        ret = {parent}
        while parent:
            parent = self.__getitem__(el)
            ret.add(parent.id)
        return ret

    def build_child_ancestor_map(self):
        ret = {tax.id: [] for tax in self.taxes}
        self.taxes.sort()
        return {tax.id: tuple(map(operator.get_attr('id'), get_parents(tax)))
                for tax in self.taxes}

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


class TaxBitMap(object):

    def __init__(self, ntaxids, k=31):
        self.map = collections.defaultdict(int)
        self.taxentries = []
        self._ntaxids = ntaxids
        self.k = k

    def add_genome(self, taxentry, path):
        assert isinstance(taxentry, TaxEntry)
        assert len(TaxEntry.genome_paths) >= 1
        self.taxentries.append(taxentries)
        assert len(self.taxentries) <= self._ntaxids
        shifted_bit = (1 << (len(self.taxentries) - 1))
        cfi = itertools.chain.from_iterable
        for kmer in set(cfi(kmer.genome2kmergen(path, k) for
                            path in taxentry.genome_paths)):
            if kmer in self.map:
                self.map[kmer] |= shifted_bit
            else:
                self.map[kmer] = shifted_bit
        sys.stderr.write(
            "Number of taxentries: %i. Bitcount of kmer %i: %i\n" % (
                len(self.taxentries), self.map[kmer],
                bin(self.map[kmer]).count("1")))

    def __str__(self):
        return ", ".join(map(str, [
            self.map, "|".join(map(str, self.taxentries)),
            self._ntaxids, self.k]))


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
                el2 = eval("%s" % el)
                assert el == el2
            els = set(tax)
            if 0 not in tax:
                print("Failed by not having 0 in tax!")
            self.assertTrue(0 in tax)
            tax.build_child_ancestor_map()
            for i in range(len(tax) - 1):
                assert tax.taxes[i].depth <= tax.taxes[i + 1]
                try:
                    assert tax.taxes[i].parent in tax
                except AssertionError:
                    print("Parent %i for %s not in tax." %
                          (tax.taxes[i].parent, tax.taxes[i]))
                    raise

    unittest.main()


__all__ = [TaxEntry, Taxonomy, load_nameidmap]
