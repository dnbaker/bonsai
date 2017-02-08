#!/usr/bin/env python
import sys
import multiprocessing
import gzip
import os
import argparse
from subprocess import check_call as cc, CalledProcessError
from enum import IntEnum
from sys import argv, stderr

class ExitCodes(IntEnum):
    EXIT_SUCCESS = 0
    EXIT_FAILURE = 1

ARCHIVE_BASE = "ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/"

GI2TAX_MAP_PATH = "ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_nucl.dmp.gz"

FOLDER = "old_ref"


def get_gi2tax(folder="oldref"):
    if not os.path.isdir(folder):
        os.makedirs(folder)
    target = folder + "/gi_taxid_nucl.dmp"
    if not os.path.isfile(target):
        cstr = "curl " + GI2TAX_MAP_PATH + " | gzip -dc > " + target
        print(cstr)
        cc(cstr, shell=True)
    return target


def parse_gi2tax(path):
    ret = {}
    for line in open(path):
        toks = line.split()
        ret[int(toks[0])] = int(toks[1])
    return ret


def getopts():
    a = argparse.ArgumentParser()
    a.add_argument("--folder", default=FOLDER)
    return a.parse_args()


def main():
    args = getopts()
    gi2tax = get_gi2tax(args.folder)
    gi2tax_map = parse_gi2tax(gi2tax)
    
    
if __name__ == "__main__":
    sys.exit(main())
