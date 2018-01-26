import sys
import shlex
import subprocess
import multiprocessing
import os


def submit_call(tup):
    k, n, opath, paths, redo = tup
    if any(not os.path.isfile(path) for path in paths):
        raise Exception("The files are NOT in the computer! %s" %
                        ", ".join(paths))
    if not redo:
        if os.path.isfile(opath) and os.path.getsize(opath):
            print("%s has run and redo is set to false. "
                  "Continuing" % opath, file=sys.stderr)
            return
    cstr = "distcmp -o%s -mn%i -k%i %s" % (opath, n, k, ' '.join(paths))
    print("Calling '%s'" % cstr, file=sys.stderr)
    subprocess.check_call(shlex.split(cstr))


def makefn(x, y, z):
    from functools import reduce
    from operator import xor
    return "experiment_%i_%x_genomes.k%i.n%i.out" % (len(x), reduce(xor, map(hash, x)),  y, z)


def main():
    sketch_range = range(14, 24, 2)
    kmer_range = range(27, 33, 2)
    import argparse
    argv = sys.argv
    # Handle args
    p = argparse.ArgumentParser(
        description="This calculates all pairwise distances between "
                    "genomes for %i combinations of parameters."
                    "This does take a while." % (len(sketch_range) *
                                                 len(kmer_range)))
    p.add_argument("--no-redo", "-n", action="store_true")
    p.add_argument("--threads", "-p",
                   default=multiprocessing.cpu_count())
    p.add_argument('genomes', metavar='paths', type=str, nargs='+',
                   help=('paths to genomes or a path to a file'
                         ' with one genome per line.'))
    args = p.parse_args()

    threads = args.threads
    redo = not args.no_redo
    paths = genomes = args.genomes
    if len(genomes) == 1 and os.path.isfile(next(open(genomes[0])).strip()):
        paths = [i.strip() for i in open(genomes[0])]
    submission_sets = ((ks, ss, makefn(paths, ks, ss), paths, redo)
                       for ss in sketch_range for ks in kmer_range)
    multiprocessing.Pool(multiprocessing.cpu_count()).map(submit_call,
                                                          submission_sets)
    return 0


if __name__ == "__main__":
    sys.exit(main())
