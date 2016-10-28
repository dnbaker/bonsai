#!/usr/bin/env python
import sys
import multiprocessing
import gzip
from subprocess import check_call as cc, CalledProcessError
argv = sys.argv


def isvalid_gzip(fn):
    try:
        cc("gunzip -t " + fn, shell=True)
        return True
    except CalledProcessError:
        print("Corrupted file ", fn, ". Delete, try again.")
        return False


def xfirstline(fn):
    # Works on python3, not 2.
    if not isvalid_gzip(fn):
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



ALL_CLADES_MAP = {
    "archaea": "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/",
    "bacteria": "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/",
    "fungi": "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/",
    "viral": "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/",
    "plant": "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/plant/",
    "protozoa": "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/protozoa/",
    "human": "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens",
    "vertebrate_mammalian": "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/",
    "vertebrate_other": "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_other/"
}

DEFAULT_CLADES = [
    "archaea", "bacteria", "viral", "human"
]

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


def get_assembly_txt(url):
    return url + "/assembly_summary.txt"


def parse_assembly(fn, fnidmap):
    to_fetch = []
    for line in open(fn):
        if line[0] == '#':
            continue
        s = line.split("\t")
        if len(s) < 14:
            print(s)
            raise Exception("Not long enough")
        if ("latest" not in line
                or (("Complete Genome" not in line and
                                     "GRCh" not in line
                     and s[13] != "Full"))):
            print("Failing line %s" % line[:-1])
            continue
        fn = "%s_genomic.fna.gz" % ([i for i in s[-2].split("/") if i][-1])
        fnidmap[fn] = int(s[5])
        to_fetch.append(s[-2] + "/" + fn)
    return to_fetch


def retry_cc(cstr):
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
                raise Exception("Could not download via %s "
                                "even after %i attempts." % (cstr,
                                                             RETRY_LIMIT))
            continue
    print("Success with %s" % cstr)


def getopts():
    import argparse
    a = argparse.ArgumentParser()
    a.add_argument("--idmap", "-m", help="Path to nameidmap.",
                   default="nameidmap.txt")
    a.add_argument("--ref", "-r", help="Name of folder for references.")
    a.add_argument("clades", nargs="+", help="Clades to use.")
    a.add_argument("--threads", "-p",
                   help="Number of threads to use while downloading.",
                   type=int, default=16)
    return a.parse_args()

def main():
    global TAX_PATH
    tax_path = TAX_PATH  # Make global variable local
    args = getopts()
    ref = args.ref if args.ref else "ref"
    import os
    if not os.path.isdir(ref):
        os.makedirs(ref)
    clades = args.clades if args.clades else DEFAULT_CLADES
    for clade in clades:
        print(clade)
        assert clade in ALL_CLADES_MAP or clade in ["all", "default"]
    to_dl = get_clade_map(clades)
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
        for s in to_dl[clade]:
            fn = "%s/%s/%s" % (ref, clade, s.split("/")[-1])
            if os.path.isfile(fn):
                if not isvalid_gzip(fn):
                    cc("rm " + fn, shell=True)
        cstrs = [("curl %s -o %s/%s/%s" %
                 (s, ref, clade, s.split("/")[-1])) for s in to_dl[clade]
                  if not os.path.isfile("%s/%s/%s" % (ref, clade,
                                                       s.split("/")[-1]))]
        # If nodes.dmp hasn't been downloaded, grab it.
        if not os.path.isfile("%s/nodes.dmp" % ref):
            cstrs.append("curl {tax_path} -o {ref}/"
                         "taxdump.tgz && tar -zxvf {ref}/taxdump.tgz"
                         " && mv nodes.dmp {ref}/nodes.dmp".format(**locals()))
        spoool = multiprocessing.Pool(args.threads)
        spoool.map(retry_cc, cstrs)
        # Replace pathnames with seqids
        for fn in list(cladeidmap.keys()):
            path = "/".join([ref, clade, fn])
            cladeidmap[xfirstline(path).decode().split()[0][1:]] = cladeidmap[fn]
            del cladeidmap[fn]
        nameidmap.update(cladeidmap)
    print("Done with all clades")
    with open(ref + "/nameidmap.txt", "w") as f:
        fw = f.write
        for k, v in nameidmap.items():
            fw(k + "\t" + str(v) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
