#!/usr/bin/env python
import sys
import multiprocessing
from subprocess import check_call as cc, CalledProcessError
argv = sys.argv



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


def get_clade_map(clades):
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


def parse_assembly(fn):
    to_fetch = []
    for line in open(fn):
        if line[0] == "#":
            continue
        s = line.split("\t")
        if len(s) < 14:
            print(s)
            raise Exception("Not long enough")
        if s[10] != "latest" or s[11] != "Complete Genome" and s[13] != "Full":
            continue
        to_fetch.append(s[-2])
    to_fetch = ["%s/%s_genomic.fna.gz" % (ftp, [i for i in ftp.split("/") if i][-1])
                for ftp in to_fetch]
    return to_fetch

def retry_cc(cstr):
    print("Starting retry_cc")
    RETRY_LIMIT = 10
    r = 0
    while r < RETRY_LIMIT:
        try:
            cc(cstr, shell=True)
        except CalledProcessError:
            r += 1
            continue
    if r == RETRY_LIMIT:
        raise Exception("Could not download via %s "
                        "even after %i attempts." % (cstr,
                                                     RETRY_LIMIT))
    print("Success with %s" % cstr)

def main():
    import os
    if not os.path.isdir("ref"):
        os.makedirs("ref")
    clades = argv[1:] if argv[1:] else DEFAULT_CLADES
    to_dl = get_clade_map(argv[1:] if argv[1:] else DEFAULT_CLADES)
    for clade in to_dl:
        if not os.path.isdir("ref/" + clade):
            os.makedirs("ref/" + clade)
        if not os.path.isfile("ref/%s/as.%s.txt" % (clade, clade)):
            cstr = ("curl %s/assembly_summary.txt "
                    "-o ref/%s/as.%s.txt") % (to_dl[clade], clade, clade)
            print(cstr)
            cc(cstr, shell=True)
        to_dl[clade] = parse_assembly("ref/%s/as.%s.txt" % (clade, clade))
        cstrs = [("curl %s -o ref/%s/%s" %
                 (assum, clade, assum.split("/")[-1])) for assum in to_dl[clade]
                  if not os.path.isfile("ref/%s/%s" % (clade,
                                                       assum.split("/")[-1]))]
        print(cstrs)
        [retry_cc(cstr) for cstr in cstrs]
    return 0


if __name__ == "__main__":
    sys.exit(main())
