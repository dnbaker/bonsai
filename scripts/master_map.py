#!/usr/bin/env python
import sys
import multiprocessing
import gzip
import os
from subprocess import check_call as cc, CalledProcessError
from download_genomes import is_valid_gzip, xfirstline
argv = sys.argv


def getopts():
    import argparse
    a = argparse.ArgumentParser()
    a.add_argument("paths", nargs="+", help="Paths to as files.")
    a.add_argument("--threads", "-p",
                   help="Number of threads to use while downloading.",
                   type=int, default=1)
    a.add_argument("-o", "--out",
                   help="Path to write output. Default: stdout")
    return a.parse_args()


def as2dict(path):
    ret = {}
    folder = "/".join(path.split("/")[:-1]) + "/"
    for line in open(path):
        if line.startswith("#"):
            continue
        toks = line.split("\t")
        taxid = int(toks[5])
        fn = folder + toks[19].split("/")[-1] + "_genomic.fna.gz"
        if not os.path.isfile(fn):
            sys.stderr.write("%s not a file. Continuing.\n" % fn)
            continue
        ret[xfirstline(fn)[1:].split()[0].decode()] = taxid
    return ret


FTP_BASENAME = "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/"


def main():
    args = getopts()
    master = {}
    mini_dicts = (multiprocessing.Pool(args.threads).map if args.threads > 1
                  else map)(as2dict, args.paths)
    with open(args.out, "w") if args.out else sys.stdout as outfile:
        of = outfile.write
        for mini in mini_dicts:
            for k, v in mini.items():
                of("%s\t%i\n" % (k, v))
    return 0


if __name__ == "__main__":
    sys.exit(main())
