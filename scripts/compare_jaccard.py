#/usr/bin/env python
import sys
import os
import multiprocessing
import subprocess
import shlex
import numpy as np

def parse_sizefile(fn):
    return {line.strip().split()[0]: float(line.strip().split()[1])
            for line in open(fn) if line[0] != "#"}

def frac_off(est, corr):
    return float(abs(est - corr)) / corr;


class DistData:
    def __init__(self, arr, names):
        self.arr = arr
        self.names = names


def parse_distfile(fn):
    fp = open(fn)
    header = next(fp)
    names = header.strip().split()[1:]
    arr = np.array([list(map(float,
                             line.strip().split()[1:]))
                    for line in fp], dtype=np.double)
    return DistData(arr, names)


def main():
    nthreads = multiprocessing.cpu_count()
    sketch_sizes = "sketch_sizes.txt"
    sketch_dists = "sketch_dists.txt"
    set_sizes = "set_sizes.txt"
    set_dists = "set_dists.txt"
    argc, argv = len(sys.argv), sys.argv
    for arg in argv:
        if "-p" in arg:
            nthreads = int(arg[2:])
            argv = [arg2 for arg2 in argv if arg2 != arg]
            print(argv)
            break
    ofp = sys.stdout
    ofw = ofp.write
    paths = argv[1:]
    for path in paths:
        if not os.path.isfile(path):
            raise RuntimeError("Path %s is not a file. Abort!" % path)
    cstr = ("bonsai dist -p%i -o %s "
            "-O %s %s") % (nthreads, sketch_sizes,
                           sketch_dists, " ".join(paths))
    print(cstr)
    subprocess.check_call(shlex.split(cstr))
    cstr = ("bonsai setdist -p%i -o %s "
            "-O %s %s") % (nthreads, set_sizes, set_dists, " ".join(paths))
    print(cstr)
    subprocess.check_call(shlex.split(cstr))
    estim_sizes = parse_sizefile(sketch_sizes)
    exact_sizes = parse_sizefile(set_sizes)
    ofw("#Name\tEstim\tExact\tError\n")
    for pair in zip(estim_sizes.items(), exact_sizes.items()):
        ofw("%s.\t%f.\t%f\t%f%%\n" % (pair[0][0], pair[0][1],
                                    pair[1][1],
                                    frac_off(pair[0][1],
                                             pair[1][1]) * 100))
    estim_dists, exact_dists = list(map(parse_distfile,
                                        [sketch_dists, set_dists]))
    names = estim_dists.names
    ofw("#Name1\tName2\tEstim\tExact\tError\n")
    for i in range(len(estim_dists.arr)):
        for j in range(i + 1, len(estim_dists.arr)):
            ofw("%s\t%s\t%f\t%f\t%f%%\n" % (names[i], names[j],
                estim_dists.arr[i,j], exact_dists.arr[i,j],
                frac_off(estim_dists.arr[i,j], exact_dists.arr[i,j]) * 100))

if __name__ == "__main__":
    sys.exit(main())
