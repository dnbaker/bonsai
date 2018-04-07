import argparse
import os
import shlex
import subprocess
import multiprocessing
import sys


def sketch_call(tup):
    ss, ks, paths = tup
    cstr = "flashdans sketch -p1 -k%i -S%i %s" % (ss, ks, " ".join(paths))
    subprocess.check_call(shlex.split(cstr))


if __name__ == "__main__":
    sketch_range = range(10, 24, 1)
    p = argparse.ArgumentParser(
        description="This calculates all pairwise distances between "
                    "genomes for all  combinations of parameters."
                    "This does take a while.")
    p.add_argument("--threads", "-p",
                   default=multiprocessing.cpu_count(), type=int)
    p.add_argument('genomes', metavar='paths', type=str, nargs='+',
                   help=('paths to genomes or a path to a file'
                         ' with one genome per line.'))
    p.add_argument("--range-start", default=24, type=int)
    p.add_argument("--range-end", default=32, type=int)
    args = p.parse_args()

    kmer_range = range(args.range_start, args.range_end + 1)
    threads = args.threads
    paths = args.genomes
    path = paths[0]
    if os.path.isfile(path) and os.path.isfile(next(open(path)).strip()):
        paths = list(open(path))
    Spooool = multiprocessing.Pool(threads)
    Spooool.map(sketch_call,
                ((ss, ks, paths) for ss in sketch_range for ks in kmer_range))
