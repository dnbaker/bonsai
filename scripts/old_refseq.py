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
    print("parsing gi2tax")
    ret = {}
    for linenum, line in enumerate(open(path)):
        if linenum % 500000 == 0:
            stderr.write("%i lines of 587716298 processed "
                         "into gi2tax map, now of size %i\n" %
                         (linenum, len(ret)))
        toks = line.split()
        taxid = int(toks[1])
        if taxid not in acceptable_taxids:
            continue
        ret[int(toks[0])] = taxid
    return ret


def append_old_to_new(gi2tax_map, newmap, concat_map, folder=FOLDER):
    paths = co("find %s -name '*.fna'" % folder, shell=True)
    if isinstance(paths, bytes):
        paths = paths.decode()
    paths = paths.split()
    ofh = open(concat_map, "w")
    ofw = ofh.write
    print("Length of accepted things: %i" % len(gi2tax_map))
    missing = set()
    for path in paths:
        sys.stderr.write("Processing path %s" % path)
        fl = xfirstline(path)
        ptoks = fl.split("|")
        name = ptoks[3]
        key = int(ptoks[1])
        try:
            ofw("%s\t%i\n" % (name, gi2tax_map[key]))
        except KeyError:
            missing.add(int(fl.split("|")[1]))
    print("Missing: " + str(missing))
    ofh.close()
    cc("cat %s >> %s" % (newmap, concat_map), shell=True)
    return concat_map


def getopts():
    a = argparse.ArgumentParser()
    a.add_argument("--folder", default=FOLDER)
    a.add_argument("--new-refseq-nameid-map", "-N", default="ref/nameidmap.txt")
    a.add_argument("--combined-nameid-map", "-c", required=True)
    a.add_argument("--taxonomy", "-t", required=True)
    a.add_argument("--no-download", '-D', default=False, action="store_true")
    return a.parse_args()


def build_full_taxmap(ncbi_tax_path):
    ret = {}
    for line in open(ncbi_tax_path):
        a, b = map(int, line.split("|")[:2])
        ret[a] = b
    return ret


def fill_set_from_tax(tax, taxmap):
    ret = {tax}
    get = taxmap.get(tax)
    while get and get > 1:
        ret.add(get)
        get = taxmap.get(get)
    return ret


def get_acceptable_taxids(taxmap):
    ret = set()
    path = "ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Bacteria/summary.txt"
    cc("curl %s > tax_summary.txt" % path, shell=True)
    for line in open("tax_summary.txt"):
        toks = line.split()
        if toks[0] == "Accession":
            continue
        tax = int(toks[3])
        if tax in ret:
            continue
        else:
            ret |= fill_set_from_tax(tax, taxmap)
            print("Size of ret: %i" % len(ret))
    print("Acceptable:", ret)
    # cc("rm tax_summary.txt", shell=True)
    return ret


def fetch_i100(folder):
    if not os.path.isdir(folder + "/i100"):
        os.makedirs(folder + "/i100")
    cstr = ("wget -N -m -np -nd -e robots=off -P %s/i100 -A .gz "
            "http://www.bork.embl.de/~mende/simulated_data/") % folder
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
    if args.no_download is False:
        fetch_genomes(args.folder)
    print("Getting acceptable taxids")
    taxmap = build_full_taxmap(args.taxonomy)
    acceptable_taxids = get_acceptable_taxids(taxmap)
    print("Appending old to new")
    concat = append_old_to_new(parse_gi2tax(gi2tax, acceptable_taxids),
                               args.new_refseq_nameid_map,
                               args.combined_nameid_map, args.folder)
    nl = int(co("wc -l %s" % concat, shell=True).decode().split()[0])
    sys.stderr.write("Concatenated file of total lines "
                     "%i is written to %s.\n" % (nl, concat))
    return 0
    
if __name__ == "__main__":
    sys.exit(main())
