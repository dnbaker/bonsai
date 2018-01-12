import sys
import shlex
import subprocess
import multiprocessing
import os


def submit_call(tup):
    k, n, opath, paths, redo = tup
    if not redo:
        if os.path.isfile(opath):
            if os.path.getsize(opath):
                print("%s has run and redo is set to false. "
                      "Continuing" % opath)
                return
    cstr = "distcmp -o%s -mn%i -k%i %s" % (opath, n, k, ' '.join(paths))
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
    if any(flag in argv for flag in flags):
        argv = [i for i in argv if i not in flags]
        redo = False
    if "-h" in argv:
        print("This calculates all pairwise distances between "
              "genomes for %i combinations of parameters."
              "This does take a while." % sum(1 for i in range(14, 24) for
                                              j in range(14, 32)))
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
                       for ss in range(14, 24) for ks in range(14, 32))
    multiprocessing.Pool(multiprocessing.cpu_count()).map(submit_call,
                                                          submission_sets)
    return 0


if __name__ == "__main__":
    sys.exit(main())
