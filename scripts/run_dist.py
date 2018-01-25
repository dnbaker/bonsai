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
    return "experiment_%i_genomes.k%i.n%i.out" % (len(x), y, z)


def usage():
    print("python %s [path to file with"
          " one fn per line] <or paths here>" % sys.argv[0])


def main():
    argv = sys.argv
    # Handle args
    flags = ("--no-redo", "-n")
    redo = True
    sketch_range = range(14, 24, 2)
    kmer_range = range(18, 32, 2)
    if any(flag in argv for flag in flags):
        argv = [i for i in argv if i not in flags]
        redo = False
    if "-h" in argv:
        print("This calculates all pairwise distances between "
              "genomes for %i combinations of parameters."
              "This does take a while." % (len(sketch_range) *
                                           len(kmer_range)),
              file=sys.stderr)
        usage()
        return 1
    try:
        paths = [i.strip() for i in open(argv[1])]
        if not all(os.path.isfile(path) for path in paths):
            raise Exception("Punt this to opening all files")
    except:
        paths = sys.argv[1:]
    if not paths:
        usage()
        return 1
    submission_sets = ((ks, ss, makefn(paths, ks, ss), paths, redo)
                       for ss in sketch_range for ks in kmer_range)
    multiprocessing.Pool(multiprocessing.cpu_count()).map(submit_call,
                                                          submission_sets)
    return 0


if __name__ == "__main__":
    sys.exit(main())
