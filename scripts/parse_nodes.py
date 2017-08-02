import sys


DEFAULT_CLADES = ["superkingdom", "kingdom", "subkingdom", "superphylum",
                  "phylum", "subphylum", "superclass", "class", "subclass",
                  "infraclass", "cohort", "superorder", "order", "suborder",
                  "infraorder", "parvorder", "superfamily", "family",
                  "subfamily", "tribe", "subtribe", "genus", "subgenus",
                  "species group", "species subgroup", "species",
                   "subspecies", "varietas", "forma"]


lmangle = lambda s: s.lower().replace(" ", "_").replace("\t", " ")
umangle = lambda s: s.upper().replace(" ", "_").replace("\t", " ")
manglers = set(lmangle, umangle)


def get_level(line):
    gen = (i for i in line.strip().split("\t") if i != "|")
    next(gen)
    next(gen)
    return next(gen)


def get_levels(path):
    from collections import Counter
    return Counter(get_level(line) for line in open(path))


def generate_enum(clades, maxl):
    enum_str = "enum class ClassLevel:int {\n"
    enum_str += "\n".join("    %s%s= %i," % (umangle(clade),
                                             (maxl - len(clade) + 1) * ' ',
                                             id)
                          for id, clade in enumerate(clades, start=2))
    enum_str += ("\n    ROOT %s= 1,"
                 "\n    NO_RANK %s= 0\n};\n\n" % ((maxl - 4) * ' ',
                                                   (maxl - 7) * ' '))
    return enum_str


def generate_class_names(clades):
    reth = "extern const char *classlvl_arr[%i];\n" % (len(clades) + 2)
    retc = "const char *classlvl_arr[%i] {\n" % (len(clades) + 2)
    retc += "    \"no rank\",\n    \"root\",\n"
    retc += "%s\n};\n\n" % "\n".join("    \"%s\"," % mangled_clade for
                                     mangled_clade in map(lmangle, clades))
    return reth, retc


def generate_class_map(clades, maxl):
    reth = ("extern const std::unordered_map"
            "<std::string, ClassLevel> classlvl_map;\n")
    retc = "const std::unordered_map<std::string, ClassLevel> classlvl_map {\n"
    retc += "    {\"no rank\", %sClassLevel::NO_RANK},\n" % ((maxl - 7) * ' ')
    retc += "    {\"root\", %sClassLevel::ROOT},\n" % ((maxl - 4) * ' ')
    for clade in clades:
        print("Adding clade level %s" % lmangle(clade))
    retc += "\n".join("    {\"%s\", %sClassLevel::%s}," %
                      (clade.lower(), (maxl - len(clade)) * ' ',
                       umangle(clade))
                      for clade in clades)
    retc += "\n};\n\n"
    return reth, retc


def generate_python_class_map(clades):
    ret = {'no rank': 0, 'root': 1}
    ret.update((el, ind) for el, ind in enumerate(
               map(lmangle, clades), start=2))
    return ret


def generate_code(clades):
    if not isinstance(clades,  list):
        try:
            clades = list(clades)
        except:
            raise RuntimeError("Clades is of type %s" % type(clades))
    clades = [clade for clade in clades if clade not in
              ("no", "rank", "no rank", "root")]
    prefix = "namespace emp {\n"
    suffix = "\n} // namespace emp\n"
    maxl = max(max(map(len, clades)), 7)
    enumstr = generate_enum(clades, maxl)
    hdr, retstr = generate_class_names(clades)
    hdr = enumstr + hdr
    chdr, cf = generate_class_map(clades, maxl)
    retstr += cf
    hdr += chdr
    retstr = prefix + retstr + suffix
    retstr = "#include \"sample_gen.h\"\n" + retstr
    hdr = "\n".join(("#ifndef _CLADE_HEADER_H__\n#define _CLADE_HEADER_H__",
                     "\n".join("#include <%s>" % i for i in
                               ("string", "iostream", "cassert",
                                "unordered_map", "cstring")),
                     prefix, hdr, "#define LINE_LVL_OFFSET 0\n", suffix,
                     "#endif // #ifndef _CLADE_HEADER_H__\n"))
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
    addstr = ""
    endstr = ""
    standalone = False
    python_outpath = ""
    for i in argv:
        if i in ["-h", "--help", "-?"]:
            sys.stderr.write("Usage: python %s <namesfile.txt> "
                             "<out.h> <out.cpp>\n"
                             "[Omit outfiles to write to "
                             "stderr/stdout respectively.]\n" % argv[0])
            sys.exit(1)
        if i in ["-c", "--standalone"]:
            endstr = ("int main() {\n    for(const auto &pair:"
                      " classlvl_map) {\n"
                      "        std::cerr << pair.first << \" has index \" << "
                      "(int)pair.second + 2 << \" and value \" << "
                      "(int)pair.second << \" and matching string \" << "
                      "classlvl_arr[(int)pair.second + 2] << '\\n';\n        "
                      "assert(std::strcmp(pair.first.data(), classlvl_arr"
                      "[(int)pair.second + 2]) == 0);}"
                      "\n}\n")
            standalone = True
        if i.split("=")[0] in ["-p", "--python"]:
           python_outpath = i.split("=")[1] if "=" in i else "a.out.py"

    argv = [i for i in argv if i not in ["-h", "-c", "--standalone",
                                         "--help", "-?"] and
            "--python=" not in i and "-p=" not in i]
    with open(argv[1]) if argv[1:] else sys.stdin as f:
        try:
            clades = [line.strip() for line in f]
        except:
            clades = DEFAULT_CLADES
    if python_outpath:
        print(generate_python_class_map(clades), file=open(python_outpath, "w"))
    doth, dotc = generate_code(clades)
    if standalone:
        addstr = "//.h:\n" + addstr
        dotc = "//.cpp:\n" + dotc
    print(addstr + doth,
          file=open(argv[2], "w") if argv[2:] else sys.stderr)
    print(dotc + endstr,
          file=open(argv[3], "w") if len(argv) >= 3 else sys.stdout)
    return 0


if __name__ == "__main__":
    sys.exit(main())
