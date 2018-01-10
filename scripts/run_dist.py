import sys
import shlex
import subprocess
import multiprocessing

def submit_call(tup):
    k, n, opath, paths = tup
    cstr = "distcmp -o%s -mn%i -k%i %s" % (opath, n, k, ' '.join(paths))
    subprocess.check_call(shlex.split(cstr))

def main():
    try:
        paths = [i.strip() for i in open(sys.argv[1])]
    except:
        paths = sys.argv[1:]
    makefn = lambda x, y, z: "experiment_%i_genomes.k%i.n%i.out" % (len(x), y, z)
    submission_sets = ((ks, ss, makefn(paths, ks, ss), paths) for ss in range(10, 28) for ks in range(14, 32, 2))
    multiprocessing.Pool(multiprocessing.cpu_count()).map(submit_call, submission_sets)
    return 0


if __name__ == "__main__":
    sys.exit(main())
