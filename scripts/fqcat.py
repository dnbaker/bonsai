#!/usr/bin/env python
import sys
from sys import stderr
import os
from subprocess import check_call as cc, CalledProcessError


def is_valid_gzip(fn):
    import gzip
    with gzip.open(fn) as f:
        try:
            f.readline()
            return True
        except:
            stderr.write("File at %s could not be read from.\n" % fn)
            raise


def fname_sort(a):
    return a.split("/")[-1].split("_")[2] + a.split("/")[-1].split("_")[4]


def make_output(paths, outpath):
    if not outpath:
        fld = "/".join(paths[0].split("/")[:-1])
        toks = paths[0].split("/")[-1].split("_")
        toks[2] = "LALL"
        toks[-1] = "ALL.fastq.gz"
        outpath = fld + "/" + "_".join(toks)
    # print("Outpath: %s" % outpath)
    cstr = "cat " + " ".join(paths) + " > " + outpath
    # print(cstr)
    cc(cstr, shell=True)
    if is_valid_gzip(outpath):
        [cc("rm " + path, shell=True) for path in paths]
    return outpath
    

def process_folder(path, fn1=None, fn2=None):
    import glob
    if not os.path.isdir(path):
        raise NotADirectoryError("Path %s is not a directory. Abort!" % path)
    stderr.write("getting paths from %s" % (path + "/*fastq.gz"))
    fqs = glob.glob(path + "/*fastq.gz")
    for fq in fqs:
        if not is_valid_gzip(fq):
            raise SystemError("File %s is not a valid gzip file. Abort!" % fq)
    r1s, r2s = [], []
    for fq in fqs:
        if "_R1_" in fq:
            r1s.append(fq)
        elif "_R2_" in fq:
            r2s.append(fq)
        else:
            raise ValueError("Fastq name '%s' with invalid "
                             "or missing read number." % fq)
    r1s.sort(key=fname_sort), r2s.sort(key=fname_sort)
    return make_output(r1s, fn1), make_output(r2s, fn2)


def main():
    argv = sys.argv
    if not argv[1:]:
        stderr.write("Usage: python " + argv[0] + "<dir> [dir2] [dir3] [...]")
        stderr.write("Concatenates fastq.gz files, verifies that they are "
                     "intact gzip files, then deletes the original fastqs.")
        return 1
    ffqs = list(map(process_folder, argv[1:]))
    for ffqpair in ffqs:
        stderr.write("Output files available at %s, %s" % ffqpair)
    return 0
    # Could be omitted by implicit return of None, which is coerced to 0

if __name__ == "__main__":
    sys.exit(main())
