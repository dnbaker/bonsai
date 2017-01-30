#!/usr/bin/env python
import sys
import multiprocessing
import gzip
import os
from subprocess import check_call as cc, CalledProcessError
argv = sys.argv


def is_valid_gzip(fn, lazy=False):
    '''
    We could instead use gunzip -t to check, but that actual requires
    iterating through the whole file, which is very slow. This is lazy,
    but at least it makes sure that it's a gzip file.
    '''
    if lazy:
        try:
            cc("zcat %s | head &>/dev/null" % fn, shell=True)
            return True
        except CalledProcessError:
            return False
    # else
    try:
        cc("gunzip -t " + fn, shell=True)
        print(fn + "is valid")
        return True
    except CalledProcessError:
        print("Corrupted file ", fn, ". Delete, try again.")
        return False
    """
    import gzip
    with gzip.open(fn) as f:
        try:
            f.readline()
            return True
        except:
            stderr.write("File at %s could not be read from.\n" % fn)
            raise
    """


def xfirstline(fn):
    # Works on python3, not 2.
    if not is_valid_gzip(fn):
        cc("rm " + fn, shell=True)
        sys.exit(main())
    if sys.version_info[0] == 3:
        if(open(fn, "rb").read(2) == b"\x1f\x8b"):
            return gzip.open(fn).readline()
    else:
        if(open(fn, "rb").read(2) == "\x1f\x8b"):
            return gzip.open(fn).readline()
    try:
        return next(open(fn))
    except StopIteration:
        cc("rm " + fn, shell=True)
        sys.exit(main())


FTP_BASENAME = "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/"

ALL_CLADES_MAP = {
    "archaea": FTP_BASENAME + "archaea/",
    "bacteria": FTP_BASENAME + "bacteria/",
    "fungi": FTP_BASENAME + "fungi/",
    "viral": FTP_BASENAME + "viral/",
    "plant": FTP_BASENAME + "plant/",
    "protozoa": FTP_BASENAME + "protozoa/",
    "human": FTP_BASENAME + "vertebrate_mammalian/Homo_sapiens",
    "vertebrate_mammalian": FTP_BASENAME + "vertebrate_mammalian/",
    "vertebrate_other": FTP_BASENAME + "vertebrate_other/"
}

DEFAULT_CLADES = [
    "archaea", "bacteria", "viral", "human"
]

DEFAULT_CLADES_STR = ", ".join(DEFAULT_CLADES)
ALL_CLADES_STR = ", ".join(ALL_CLADES_MAP.keys())

TAX_PATH = "ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"


def get_clade_map(clades):
    if clades[0] == "default":
        return {k: v for k, v in ALL_CLADES_MAP.items() if k in DEFAULT_CLADES}
    if clades[0] == "all":
        return ALL_CLADES_MAP
    ret = {}
    clades = [i.lower() for i in clades]
    for clade in clades:
        if clade in ALL_CLADES_MAP:
            ret[clade] = ALL_CLADES_MAP[clade]
        else:
            raise ValueError("clade %s not available. Abort!" % clade)
    return ret


def parse_assembly(fn, fnidmap):
    to_fetch = []
    for line in open(fn):
        if line[0] == '#':
            continue
        s = line.split("\t")
        if len(s) < 14:
            print(s)
            raise Exception("Not long enough")
        if ("latest" not in line or
                (("Complete Genome" not in line and
                  "GRCh" not in line and s[13] != "Full")) or
                any(i in line.lower() for i in ["supercontig", "scaffold"])):
            continue
        fn = "%s_genomic.fna.gz" % ([i for i in s[19].split("/") if i][-1])
        fnidmap[fn] = int(s[5])
        to_fetch.append(s[-2] + "/" + fn)
    return to_fetch


def retry_cc(cstr, die=True):
    print("Starting retry_cc")
    RETRY_LIMIT = 10
    r = 0
    while r < RETRY_LIMIT:
        try:
            print(cstr)
            cc(cstr, shell=True)
            return
        except CalledProcessError:
            print("retry number", r)
            r += 1
            if r == RETRY_LIMIT:
                if die:
                    raise Exception(
                        "Could not download via %s "
                        "even after %i attempts." % (cstr, RETRY_LIMIT))
                else:
                    sys.stderr.write(
                        "Could not download %s even after %i attempts" % (
                            cstr, RETRY_LIMIT))
                return
            continue
    print("Success with %s" % cstr)


def getopts():
    import argparse
    a = argparse.ArgumentParser()
    a.add_argument("--idmap", "-m", help="Path to which to write nameidmap.",
                   default="nameidmap.txt")
    a.add_argument("--ref", "-r", help="Name of folder for references.")
    a.add_argument("clades", nargs="+", help="Clades to use."
                   " default includes %s. all includes %s." % (
                        DEFAULT_CLADES_STR, ALL_CLADES_STR))
    a.add_argument("--threads", "-p",
                   help="Number of threads to use while downloading.",
                   type=int, default=16)
    a.add_argument("--lazy", "-l", default=False, type=bool,
                   help="Don't check full gzipped file contents.")
    return a.parse_args()


def check_path(fn, lazy=False):
    print("Checking path " + fn)
    if os.path.isfile(fn):
        if not is_valid_gzip(fn, lazy=lazy):
            cc("rm " + fn, shell=True)


def check_path_lazy(path):
    check_path(path, lazy=True)


def main():
    global TAX_PATH
    tax_path = TAX_PATH  # Make global variable local
    args = getopts()
    ref = args.ref if args.ref else "ref"
    if not os.path.isdir(ref):
        os.makedirs(ref)
    clades = args.clades if args.clades else DEFAULT_CLADES
    for clade in clades:
        try:
            assert clade in ALL_CLADES_MAP or clade in ["all", "default"]
        except AssertionError:
            print("Clade %s not 'all', 'default', or one of the valid "
                  "clades: %s" % (clade, ALL_CLADES_STR))
            sys.exit(1)
    to_dl = get_clade_map(clades)
    print("About to download clades %s" % ", ".join(to_dl))
    nameidmap = {}
    for clade in to_dl:
        cladeidmap = {}
        if not os.path.isdir(ref + "/" + clade):
            os.makedirs(ref + "/" + clade)
        if not os.path.isfile("%s/%s/as.%s.txt" % (ref, clade, clade)):
            cstr = ("curl %s/assembly_summary.txt "
                    "-o %s/%s/as.%s.txt") % (to_dl[clade], ref, clade, clade)
            print(cstr)
            cc(cstr, shell=True)
        to_dl[clade] = parse_assembly("%s/%s/as.%s.txt" %
                                      (ref, clade, clade), cladeidmap)
        spoool = multiprocessing.Pool(args.threads)
        spoool.map(check_path_lazy if args.lazy else check_path,
                   ("/".join([ref, clade, s.split("/")[-1]]) for
                    s in to_dl[clade]))
        print(to_dl[clade])
        cstrs = [("curl %s -o %s/%s/%s" %
                 (s, ref, clade, s.split("/")[-1])) for
                 s in to_dl[clade] if not os.path.isfile(
                     "%s/%s/%s" % (ref, clade, s.split("/")[-1]))]
        # If nodes.dmp hasn't been downloaded, grab it.
        if not os.path.isfile("%s/nodes.dmp" % ref):
            cstrs.append("curl {tax_path} -o {ref}/"
                         "taxdump.tgz && tar -zxvf {ref}/taxdump.tgz"
                         " && mv nodes.dmp {ref}/nodes.dmp".format(**locals()))
        spoool.map(retry_cc, cstrs)
        # Replace pathnames with seqids
        for fn in list(cladeidmap.keys()):
            cladeidmap[xfirstline("/".join(
                [ref, clade, fn]
            )).decode().split()[0][1:]] = cladeidmap[fn]
            del cladeidmap[fn]
        nameidmap.update(cladeidmap)
    print("Done with all clades")
    with open(ref + "/" + args.idmap, "w") as f:
        fw = f.write
        for k, v in nameidmap.items():
            fw(k + "\t" + str(v) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
