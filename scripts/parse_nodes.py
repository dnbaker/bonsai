import sys


def get_level(line):
    gen = (i for i in line.strip().split("\t") if i != "|")
    next(gen)
    next(gen)
    return next(gen)


def get_levels(path):
    from collections import Counter
    return Counter(get_level(line) for line in open(path))


def generate_enum(clades, maxl):
    enum_str = "enum class ClassLvl:int {\n"
    enum_str += "\n".join("    %s%s= %i," % (clade.upper(),
                                             (maxl - len(clade) + 1) * ' ',
                                             id)
                          for id, clade in enumerate(clades))
    enum_str += ("\n    ROOT %s= -1,"
                 "\n    NO_RANK %s= -2\n};\n\n" % ((maxl - 4) * ' ',
                                                   (maxl - 7) * ' '))
    return enum_str


def generate_class_names(clades):
    reth = "extern const char *classlvl_arr[%i];\n" % (len(clades) + 2)
    retc = "const char *classlvl_arr[%i] {\n" % (len(clades) + 2)
    retc += "    \"no rank\",\n    \"root\"\n"
    retc += "\n".join("    \"%s\"," % clade.lower() for clade in clades)
    retc += "\n};\n\n"
    return reth, retc


def generate_class_map(clades, maxl):
    reth = ("extern const std::unordered_map"
            "<std::string, ClassLevel> classlvl_map;\n")
    retc = "const std::unordered_map<std::string, ClassLevel> classlvl_map {\n"
    retc += "    {\"no rank\", %sClassLevel::NO_RANK},\n" % ((maxl - 7) * ' ')
    retc += "    {\"root\", %sClassLevel::ROOT},\n" % ((maxl - 4) * ' ')
    retc += "\n".join("    {\"%s\", %sClassLevel::%s}," %
                      (clade.lower(), (maxl - len(clade)) * ' ',
                       clade.upper())
                      for clade in clades)
    retc += "\n};\n\n"
    return reth, retc


def generate_code(clades):
    if not isinstance(clades,  list):
        try:
            clades = list(clades)
        except:
            raise RuntimeError("Clades is of type %s" % type(clades))
    clades = [clade for clade in clades if clade not in
              ("no", "rank", "no rank", "root")]
    maxl = max(max(map(len, clades)), 7)
    enumstr = generate_enum(clades, maxl)
    hdr, retstr = generate_class_names(clades)
    hdr = enumstr + hdr
    chdr, cf = generate_class_map(clades, maxl)
    retstr += cf
    hdr += chdr
    return hdr, retstr


def print_levels():
    import sys
    if len(sys.argv) < 2:
        print("Usage: python " + sys.argv[0] + " <nodes.dmp>")
        sys.exit(1)
    print(repr(get_levels(sys.argv[1])))


def main():
    import sys
    argv = sys.argv
    for i in argv:
        if i in ["-h", "--help", "-?"]:
            sys.stderr.write("Usage: python %s <namesfile.txt> "
                             "<out.h> <out.cpp>\n"
                             "[Omit outfiles to write to "
                             "stderr/stdout respectively.]\n" % argv[0])
            sys.exit(1)
    with open(argv[1]) if argv[1:] else stdout as f:
        clades = [line.strip() for line in f]
    doth, dotc = generate_code(clades)
    print("//.h:\n" + doth,
          file=open(argv[2], "w") if argv[2:] else sys.stderr)
    print("//.cpp:\n" + dotc,
          file=open(argv[3], "w") if len(argv) >= 3 else sys.stdout)
    return 0


'''
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
Sorted:
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
    sys.exit(main())
