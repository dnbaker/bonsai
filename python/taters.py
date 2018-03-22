#!/usr/bin/env python
import sys
import os
from functools import reduce  # Oh, come on, Python 3. Reduce is awesome.
if sys.version_info[0] < 3:
    raise ImportError("Cannot import %s. Python 3+ required. "
                      "[Make your life easier, take the plunge!]" %
                      os.path.basename(__file__))


def nk2sketchbytes(sketchtype, n=0, k=0):
    if not (k and n):
        raise Exception("Set n, k kwargs for nk2sketchbytes pls")
    try:
        return {"mash": 8 if k > 16 else 4,
                "hlash": 1}[sketchtype.lower()] << n
    except KeyError:
        sys.stderr.write(
            ("Key error in nk2sketchbytes."
             "key: %s. Expected mash or hlash\n" % sketchtype))
        raise


def get_flag(tok):
    '''Also works as a boolean, as these flags are guaranteed to be nonzero
    '''
    if not tok:
        return False
    flagchar, rest = tok[0], tok[1:]
    if rest:
        try:
            return int(rest)
        except ValueError:
            # sys.stderr.write("Could not parse flag '%s'\n" % tok)
            pass
    return False


def canonname(k1, k2=None):
    if not k2:
        assert len(k1) == 2
        tmp1, k2 = k1
        k1 = tmp1
    return sorted(map(os.path.basename, (k1, k2)))


def get_nk(path):
    k = n = 0
    for tok, val in ((tok, get_flag(tok)) for tok in
                     os.path.basename(path).split(".")):
        if not val:
            continue
        if tok[0] == 'k':
            k = val
        elif tok[0] == 'n':
            n = val
        elif "exper" not in tok:
            raise ValueError("Unexpected token %s" % tok)
    if not (k and n):
        raise Exception("Invalid mashexp filename %s" % path)
    return n, k


def canonkey(k1, k2=None, n=0, k=0):
    if isinstance(k1, HlashData):
        return k1.key
    if not (k and n):
        raise Exception("ZOMG")
    return tuple((*canonname(k1, k2), n, k))


def mash2dict(path):
    with open(path) as ifp:
        n, k = get_nk(path)
        ret = {}
        for line in ifp:
            toks = line.split()
            ret[canonkey(toks[0], toks[1], nk2sketchbytes("mash", n, k), k)] \
                = float(toks[2])
    return ret


def folder2paths(paths, reqstr=""):
    # This is a misnomer. It holds either a list of paths or a string
    # representing a folder from which we glob everything
    if isinstance(paths, str):
        if os.path.isdir(paths):
            import glob
            paths = glob.glob(paths + "/*out")
    assert isinstance(paths, list)
    if len(paths) == 1:
        import glob
        paths = list(glob.glob("%s/*" % paths[0]))
    fail = False
    for path in paths:
        if not os.path.isfile(path):
            print("path %s is not a file." % path, file=sys.stderr)
            fail = True
    if fail:
        raise Exception("The files are NOT in the computer.")
    if reqstr:
        paths = [i for i in paths if reqstr in i]
    return paths


def folder2mashdict(paths):
    return reduce(lambda x, y: {**x, **mash2dict(y)},
                  folder2paths(paths, "experim"), {})


def mashdict2tsv(md, path):
    fp = path if hasattr(path, 'write') else open(path, "w")
    fw = fp.write
    fw("#Path1\tPath2\tNumber of elements\tKmer size"
       "\tNumber of bytes in sketch\tEstimated Jaccard Index\n")
    set(map(lambda k: not fw("%s\t%s\t%i\t%i\t%i\t%f\n" %
                             (k[0][0], k[0][1], k[0][2], k[0][3],
                              k[0][2] * (8 if k[0][3] > 16 else 4), k[1])),
            md.items()))
    if fp != sys.stdout:
        fp.close()


def hlashdict2tsv(hd, path):
    fp = path if hasattr(path, 'write') else open(path, "w")
    fw = fp.write
    fw("#Path1\tPath2\tSketch p\tKmer size"
       "\tNumber of bytes in sketch\tEstimated Jaccard Index"
       "\tExact Jaccard Index\n")
    set(map(lambda v: not fw(str(v)), hd.values()))
    if fp != sys.stdout:
        fp.close()


class HlashData:
    def __init__(self, line):
        toks = line.strip().split()
        self.paths = canonname(toks[:2])
        self.est, self.exact = map(float, toks[2:4])
        self.n, self.k = map(int, toks[-2:])
        self.nbytes = 1 << self.n
        self.key = canonkey(self.paths[0], self.paths[1],
                            n=self.nbytes, k=self.k)

    def __str__(self):
        return "%s\t%s\t%i\t%i\t%i\t%f\t%f\n" % (*self.paths, self.n, self.k,
                                                 self.nbytes, self.est,
                                                 self.exact)


def folder2hlashdict(paths):
    paths = folder2paths(paths, "genome")
    return {canonkey(x): x for x in map(HlashData,
                                        (line for path in paths
                                         for line in open(path) if
                                         line[0] != "#"))}


if __name__ == "__main__":
    import argparse
    import sys
    p = argparse.ArgumentParser()
    p.add_argument('subcommand', choices=('mash', 'hlash'),
                   help="Whether the parsed experiment is mash or hlash.")
    p.add_argument("path", help="Path to folder to parse files from.")
    p.add_argument("outpath", help="Path to folder to write.",
                   type=argparse.FileType('w'), default=sys.stdout)
    args = p.parse_args()
    path, outpath = args.path, args.outpath
    if args.subcommand == "mash":
        mashdict2tsv(folder2mashdict(path), outpath)
    elif args.subcommand == "hlash":
        hlashdict2tsv(folder2hlashdict(path), outpath)
    else:
        raise NotImplementedError("This never happens because choice "
                                  "is an illusion.")
    sys.exit(0)


__all__ = ["mash2dict", "folder2mashdict", "get_flag",
           "folder2hlashdict", "folder2paths",
           "HlashData", "nk2sketchbytes", "canonname", "canonkey"]
