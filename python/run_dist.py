import sys
import shlex
import subprocess
import multiprocessing
import os



def is_nonempty_file(path):
    return os.path.isfile(path) and os.path.getsize(path)


def submit_distcmp_call(tup):
    k, n, opath, paths, redo = tup
    if any(not os.path.isfile(path) for path in paths):
        raise Exception("The files are NOT in the computer! %s" %
                        ", ".join(paths))
        if not redo and is_nonempty_file(opath):
            print("%s has run and redo is set to false. "
                  "Continuing" % opath, file=sys.stderr)
            return
    cstr = "distcmp -o%s -mn%i -k%i %s" % (opath, n, k, ' '.join(paths))
    print("Calling '%s'" % cstr, file=sys.stderr)
    subprocess.check_call(shlex.split(cstr))


def makefn(x, y, z, mashify):
    from functools import reduce
    from operator import xor
    hashval = reduce(xor, map(hash, x))
    if not mashify:
        return "experiment_%i_%x_genomes.k%i.n%i.out" % (len(x), hashval, y, z)
    return "mashed_experiment_%i_%x_genomes.k%i.n%i." % (
        len(x), hashval,  y, z)


def main():
    sketch_range = range(10, 24, 1)
    import argparse
    argv = sys.argv
    # Handle args
    p = argparse.ArgumentParser(
        description="This calculates all pairwise distances between "
                    "genomes for all  combinations of parameters."
                    "This does take a while.")
    p.add_argument("--no-redo", "-n", action="store_true")
    p.add_argument("--threads", "-p",
                   default=multiprocessing.cpu_count(), type=int)
    p.add_argument("--use-mash", "-M", action="store_true",
                   help=("Use Mash to calculate distances "
                         "rather than 'bonsai dist'"))
    p.add_argument('genomes', metavar='paths', type=str, nargs='+',
                   help=('paths to genomes or a path to a file'
                         ' with one genome per line.'))
    p.add_argument("--range-start", default=24, type=int)
    p.add_argument("--range-end", default=32, type=int)
    args = p.parse_args()

    kmer_range = range(args.range_start, args.range_end + 1)
    threads = args.threads
    redo = not args.no_redo
    paths = genomes = args.genomes
    mashify = args.use_mash
    if len(genomes) == 1 and os.path.isfile(next(open(genomes[0])).strip()):
        paths = [i.strip() for i in open(genomes[0])]
    submission_sets = ((ks, ss, makefn(paths, ks, ss, mashify), paths, redo)
                       for ss in sketch_range for ks in kmer_range)
    if mashify:
        for ss in sketch_range:
            for ks in kmer_range:
                fn = makefn(paths, ks, ss, mashify)
                for i in range(len(paths) - 1):
                    identifier = os.path.basename(paths[i]).split(".")[0]
                    thisfn = "%s.%s.out" % (fn, identifier)
                    if not redo and is_nonempty_file(thisfn):
                        continue
                    cstr = "mash dist -s %i -t -k %i -p %i %s > %s" % (
                        1 << ss, ks, threads,
                        " ".join(paths[i:]), thisfn)
                    subprocess.check_call(cstr, shell=True)
    else:
        multiprocessing.Pool(threads).map(submit_distcmp_call,
                                          submission_sets)
    return 0


if __name__ == "__main__":
    sys.exit(main())
