#!/usr/bin/env python3
import sys
import multiprocessing
import gzip
import os
from subprocess import check_call as cc, CalledProcessError
from enum import IntEnum
argv = sys.argv


if sys.version_info[0] != 3:
    raise Exception("Python 3 required")


class ExitCodes(IntEnum):
    EXIT_SUCCESS = 0
    EXIT_FAILURE = 1


def is_valid_gzip(fn, lazy=False, use_pigz=False):
    '''
    We could instead use gunzip -t to check, but that actual requires
    iterating through the whole file, which is very slow. This is lazy,
    but at least it makes sure that it's a gzip file.

    lazy simply tries to see if the first 10 lines can be read.
    It isn't very safe.

    use_pigz uses pigz instead of gzip. A bad idea if a number of processes
    have already been spawned.
    '''
    if lazy:
        try:
            cc("gzip -dc %s | head &>/dev/null" % fn, shell=True)
            return True
        except CalledProcessError:
            return False
    # lazy has already returned. This is the "else".
    cmd = ("pigz" if use_pigz else "gzip") + " -dc "
    try:
        cc(cmd + " -t " + fn, shell=True)
        sys.stderr.write(fn + " is valid\n")
        return True
    except CalledProcessError:
        sys.stderr.write("Corrupted file " + fn + ". Delete, try again.\n")
        return False


def xfirstline(fn):
    # Works on python3, not 2.
    ffn = gzip.open if open(fn, "rb").read(2) == b"\x1f\x8b" else open
    with ffn(fn) as f:
        ret = next(f)
    return ret


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
    "vertebrate_other": FTP_BASENAME + "vertebrate_other/",
    "invertebrate": FTP_BASENAME + "invertebrate/"
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
    # print(fn)
    to_fetch = []
    for line in open(fn, encoding='utf8'):
        if line[0] == '#':
            continue
        s = line.split("\t")
        if len(s) < 14:
            print(s)
            raise Exception("Not long enough")
        if ("latest" not in line or  # Complete genome
                (("Complete Genome" not in line and
                  "GRCh" not in line and "Full" not in line)) or
                any(i in line.lower() for
                    i in ["contig", "supercontig"])):
            continue
        #print(s[19], file=sys.stderr)
        fn = "%s_genomic.fna.gz" % ([i for i in s[19].split("/") if i][-1])
        fnidmap[fn] = int(s[5])
        index = len(s) - 1
        while "ftp" not in s[index] and index > 0:
            index = index - 1
        if index:
            to_fetch.append(s[index] + "/" + fn)
        else:
            print(f"No link found for {fn}, continue", file=sys.stderr)
            continue
            #raise RuntimeError("ftp link not found. line: %s" % line[:-1])
    return to_fetch


def retry_cc(tup):
    cstr, die = tup
    RETRY_LIMIT = 10
    r = 0
    while r < RETRY_LIMIT:
        try:
            print(cstr, file=sys.stderr)
            cc(cstr, shell=True)
            return
        except CalledProcessError:
            print("retry number", r, file=sys.stderr)
            r += 1
    if die:
        raise Exception(
            "Could not download via %s "
            "even after %i attempts." % (cstr, RETRY_LIMIT))
    else:
        sys.stderr.write(
            "Could not download %s even after %i attempts" % (
                cstr, RETRY_LIMIT))


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
    a.add_argument("--lazy", "-l", action='store_true',
                   help="Don't check full gzipped file contents.")
    a.add_argument("--die", "-d", action='store_true')
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
    if argv[1:] and argv[1] == "nodes":
        if not os.path.isfile("%s/nodes.dmp" % ref):
            cc("curl -s {tax_path} -o {ref}/"
               "taxdump.tgz && tar -zxvf {ref}/taxdump.tgz"
               " && mv nodes.dmp {ref}/nodes.dmp".format(**locals()),
                shell=True)
            return 0
    if not os.path.isdir(ref):
        os.makedirs(ref)
    clades = args.clades if args.clades else DEFAULT_CLADES
    for clade in clades:
        try:
            assert clade in ALL_CLADES_MAP or clade in ["all", "default"]
        except AssertionError:
            print("Clade %s not 'all', 'default', or one of the valid "
                  "clades: %s" % (clade, ALL_CLADES_STR), file=sys.stderr)
            sys.exit(ExitCodes.EXIT_FAILURE)
    to_dl = get_clade_map(clades)
    print("About to download clades %s" % ", ".join(to_dl), file=sys.stderr)
    nameidmap = {}
    for clade in to_dl:
        cladeidmap = {}
        if not os.path.isdir(ref + "/" + clade):
            os.makedirs(ref + "/" + clade)
        if not os.path.isfile("%s/%s/as.%s.txt" % (ref, clade, clade)):
            cstr = ("curl -s %s/assembly_summary.txt "
                    "-o %s/%s/as.%s.txt") % (to_dl[clade], ref, clade, clade)
            # print(cstr)
            cc(cstr, shell=True)
        to_dl[clade] = parse_assembly("%s/%s/as.%s.txt" %
                                      (ref, clade, clade), cladeidmap)
        print("Performing parallel download of clade " + str(clade), file=sys.stderr, flush=True)
        with multiprocessing.Pool(args.threads) as spoool:
            spoool.map(check_path_lazy if args.lazy else check_path,
                       ("/".join([ref, clade, s.split("/")[-1]]) for
                        s in to_dl[clade]))
        print("Checked existing paths for download of clade " + str(clade), file=sys.stderr, flush=True)
        cstrs = [("curl -s %s -o %s/%s/%s" %
                 (s, ref, clade, s.split("/")[-1])) for
                 s in to_dl[clade] if not os.path.isfile(
                     "%s/%s/%s" % (ref, clade, s.split("/")[-1]))]
        # If nodes.dmp hasn't been downloaded, grab it.
        if not os.path.isfile("%s/nodes.dmp" % ref):
            cstrs.append("curl -s {tax_path} -o {ref}/"
                         "taxdump.tgz && tar -zxvf {ref}/taxdump.tgz"
                         " && mv nodes.dmp {ref}/nodes.dmp".format(**locals()))
        print("Performing download of clade " + str(clade), file=sys.stderr, flush=True)
        with multiprocessing.Pool(args.threads) as spoool:
            spoool.map(retry_cc, ((cs, args.die) for cs in cstrs))
        print("Performed download of clade " + str(clade), file=sys.stderr, flush=True)
        # Replace pathnames with seqids
        print("Updating cladeidmap", file=sys.stderr, flush=True)
        nk = len(cladeidmap.keys())
        for i, fn in enumerate(list(cladeidmap.keys())):
            try:
                #print(ref, clade, fn)
                cladeidmap[xfirstline("/".join(
                    [ref, clade, fn]
                )).decode().split()[0][1:]] = cladeidmap[fn]
                del cladeidmap[fn]
            except FileNotFoundError:
                if args.die:
                    raise
            if (i & 0xFFF) == 0:
                print(f"{i}/{nk} in clade processed.\n", file=sys.stderr, flush=True)
        nameidmap.update(cladeidmap)
        print("Finished clade " + clade, file=sys.stderr, flush=True)
    print("Done with all clades", file=sys.stderr, flush=True)
    with open(ref + "/" + args.idmap, "w") as f:
        fw = f.write
        for k, v in nameidmap.items():
            fw(k + "\t" + str(v) + "\n")
    return ExitCodes.EXIT_SUCCESS


if __name__ == "__main__":
    sys.exit(main())
