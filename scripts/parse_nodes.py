import sys


def get_level(line):
    gen = (i for i in line.strip().split("\t") if i != "|")
    next(gen)
    next(gen)
    return next(gen)


def get_levels(path):
    from collections import Counter
    return Counter(get_level(line) for line in open(path))



'''
  --- Changed:
         'family': 8578,
         'forma': 470,
         'genus': 82587,
         'no rank': 214584,
         'order': 1418,
         'subspecies': 20777,
         'varietas': 6789


Counter({'class': 309,
         'cohort': 3,
         'family': 8575,
         'forma': 469,
         'genus': 82531,
         'infraclass': 15,
         'infraorder': 104,
         'kingdom': 3,
         'no rank': 214552,
         'order': 1421,
         'parvorder': 10,
         'phylum': 219,
         'species': 1266136,
         'species group': 438,
         'species subgroup': 133,
         'subclass': 132,
         'subfamily': 2763,
         'subgenus': 1320,
         'subkingdom': 1,
         'suborder': 331,
         'subphylum': 28,
         'subspecies': 20773,
         'subtribe': 479,
         'superclass': 4,
         'superfamily': 805,
         'superkingdom': 5,
         'superorder': 49,
         'superphylum': 2,
         'tribe': 1996,
         'varietas': 6781})
superkingdom
kingdom
subkingdom
superphylum
phylum
subphylum
superclass
class
subclass
infraclass
cohort
superorder
order
suborder
infraorder
parvorder
superfamily
family
subfamily
tribe
subtribe
genus
subgenus
species group
species subgroup
species
subspecies
varietas
forma
'''


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: python " + sys.argv[0] + " <nodes.dmp>")
        sys.exit(1)
    print(repr(get_levels(sys.argv[1])))
