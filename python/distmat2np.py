import os
import numpy as np


def fsize(path):
    return os.stat(path).st_size


def parse_file(f1, f2):
    if not all(map(os.path.isfile, (f1, f2))):
        raise FileNotFoundError("The files are NOT in the computer.")
    if fsize(f1) < fsize(f2):
        f3 = f1; f1 = f2; f2 = f3
    labels = [line.strip() for line in open(f2) if line[0] != "#"]
    nl = len(labels)
    data = np.ones((nl, nl), dtype=np.float32) * 137
    ifp = open(f1)
    next(ifp)  # discard unused
    i = 0
    for line in ifp:
        toks = line.strip().split('\t')[1:]
        j = i + 1
        data[i,i] = 1.  # bc the ji of a set and itself is 1
        while j < nl:
            data[j,i] = data[i,j] = float(toks[j])
            j += 1
        i += 1
    return data


def usage():
    import sys
    sys.stderr.write(f"{sys.argv[0]} <distances> <sizes> [outfile] (default: extracted.npy)\n")
    sys.exit(1)


if __name__ == "__main__":
    import sys
    argv = sys.argv
    if len(argv) < 3:
        usage()
    data = parse_file(argv[1], argv[2])
    np.save(argv[3] if argv[3:] else "extracted.npy", data)


__all__ = [parse_file, fsize, usage]
