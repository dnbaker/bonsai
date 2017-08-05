#!/usr/bin/env python
import sys
import os
import multiprocessing
from sys import stderr
from subprocess import check_call as cc, CalledProcessError


def is_valid_gzip(fn):
    '''
    We could instead use gunzip -t to check, but that actual requires
    iterating through the whole file, which is very slow. This is lazy,
    but at least it makes sure that it's a gzip file.
    '''
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
        if not paths:
            print("Missing input in folder")
            return
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
        [cc("rm -f " + path, shell=True) for path in paths]
    return outpath


def process_folder(path, fn1=None, fn2=None):
    import glob
    if not os.path.isdir(path):
        raise NotADirectoryError("Path %s is not a directory. Abort!" % path)
    stderr.write("getting paths from %s" % (path + "/*fastq.gz"))
    fqs = [i for i in glob.glob(path + "/*fastq.gz") if "LALL" not in i]
    if not fqs:
        already_done = glob.glob(path + "/*LALL*fastq.gz")
        if already_done:
            return ([i for i in already_done if "_R1_" in i][0],
                    [i for i in already_done if "_R2_" in i][0])
        else:
            raise FileNotFoundError("Could not find any files in the folder" +
                                    path)
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
    ret = (make_output(r1s, fn1), make_output(r2s, fn2))
    if not ret[0] or not ret[1]:
        print("Could not get files from " + path)
        return None


def getopts():
    import argparse
    a = argparse.ArgumentParser(
        description="Concatenates fastq.gz files, verifies that they are "
                    "intact gzip files, then deletes the original fastqs.")
    a.add_argument("folders", nargs="+", help="Folders to process")
    a.add_argument("--threads", "-p",
                   help="Number of threads to use.",
                   type=int, default=8)
    return a.parse_args()


def main():
    opts = getopts()
    spoool = multiprocessing.Pool(opts.threads)
    ffqs = spoool.map(process_folder, opts.folders)
    for index, ffqpair in enumerate(ffqs):
        if ffqpair:
            stderr.write("Output files available at %s, %s\n" % ffqpair)
        else:
            stderr.write("Error processing directory %s\n" %
                         opts.folders[index])
    return any(not x for x in ffqs)
    # Could be omitted by implicit return of None, which is coerced to 0

if __name__ == "__main__":
    sys.exit(main())
