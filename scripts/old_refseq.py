#!/usr/bin/env python
import sys
import multiprocessing
import gzip
import os
import argparse
from subprocess import (check_call as cc, CalledProcessError,
                        check_output as co)
from enum import IntEnum
from sys import argv, stderr
from download_genomes import xfirstline

class ExitCodes(IntEnum):
    EXIT_SUCCESS = 0
    EXIT_FAILURE = 1

ARCHIVE_BASE = "ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq"

GI2TAX_MAP_PATH = "ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_nucl.dmp.gz"

FOLDER = "ref/old"

NAMES = ["Bacteria"]


def get_gi2tax(folder=FOLDER):
    if not os.path.isdir(folder):
        os.makedirs(folder)
    target = folder + "/gi_taxid_nucl.dmp"
    if not os.path.isfile(target):
        cc("curl %s | gzip -dc > %s" % (GI2TAX_MAP_PATH, target), shell=True)
    return target


def parse_gi2tax(path, acceptable_taxids=set()):
    ret = {}
    for line in open(path):
        toks = line.split()
        taxid = int(toks[1])
        if taxid not in acceptable_taxids:
            continue
        ret[int(toks[0])] = int(toks[1])
        if len(ret) % 5000000 == 0:
            stderr.write("%i members of 587716298 loaded "
                         "into gi2tax map\n" % len(ret))
    return ret


def append_old_to_new(gi2tax_map, newmap, concat_map, folder=FOLDER):
    paths = co("find %s -name '*.fna'" % folder, shell=True)
    if isinstance(paths, bytes):
        paths = paths.decode()
    paths = paths.split()
    ofh = open(concat_map, "w")
    ofw = ofh.write
    for path in paths:
        fl = xfirstline(path)
        ofw("%s\t%i\n" % (fl.split()[1], gi2tax_map[int(fl.split("|")[1])]))
    ofh.close()
    cc("cat %s >> %s" % (newmap, concat_map), shell=True)
    return concat_map


def getopts():
    a = argparse.ArgumentParser()
    a.add_argument("--folder", default=FOLDER)
    a.add_argument("--new-refseq-nameid-map", "-N", default="ref/nameidmap.txt")
    a.add_argument("--combined-nameid-map", "-c", required=True)
    return a.parse_args()


def get_acceptable_taxids():
    ret = set()
    path = "ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Bacteria/summary.txt"
    cc("curl %s > tax_summary.txt" % path, shell=True)
    for line in open("tax_summary.txt"):
        toks = line.split()
        if toks[0] == "Accession":
            continue
        ret.add(int(toks[3]))
    cc("rm tax_summary.txt", shell=True)


def fetch_i100(folder):
    if not os.path.isdir(folder + "/i100"):
        os.makedirs(folder + "/i100")
    cstr = ("wget -N -m -np -nd -e robots=off -P %s/i100 -A .gz "
            "http://www.bork.embl.de/~mende/simulated_data/")
    cc(cstr, shell=True)


def fetch_genomes(folder, names=NAMES):
    for name in names:
        ftp_path = "/".join([ARCHIVE_BASE, name])
        cstr = ("wget -N -m -np -nd -e robots=off -P"
                " %s/%s -A .fna,.fna.gz %s") % (folder, name, ftp_path)
        cc(cstr, shell=True)
    fetch_i100(folder)


def main():
    args = getopts()
    gi2tax = get_gi2tax(args.folder)
    fetch_genomes(args.folder)
    acceptable_taxids = get_acceptable_taxids()
    concat = append_old_to_new(parse_gi2tax(gi2tax), args.new_refseq_nameid_map,
                               args.combined_nameid_map, args.folder)
    nl = int(co("wc -l %s" % concat, shell=True).decode().split()[0])
    sys.stderr.write("Concatenated file of total lines "
                     "%i is written to %s.\n" % (nl, concat))
    return 0
    
if __name__ == "__main__":
    sys.exit(main())
