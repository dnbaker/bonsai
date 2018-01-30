#!/usr/bin/env python
import sys
import os
if sys.version_info[0] < 3:
    raise ImportError("Cannot import %s. Python 3+ required. "
                      "[Make your life easier, take the plunge!]" %
                      os.path.basename(__file__))

basename = os.path.basename


def get_flag(tok):
    '''Also works as a boolean, as these flags are guaranteed to be nonzero
    '''
    if not tok:
        return False
    flagchar = tok[0]
    rest = tok[1:]
    if rest:
        try:
            return int(rest)
        except ValueError:
            # sys.stderr.write("Could not parse flag '%s'\n" % tok)
            pass
    return False


def mash2dict(path):
    with open(path) as ifp:
        flags = [(tok, get_flag(tok)) for
                 tok in path.split("/")[-1].split(".")]
        k = 0
        n = 0
        print(flags)
        for tok, val in flags:
            if not val:
                continue
            if tok[0] == 'k':
                k = val
            elif tok[0] == 'n':
                n = val
            elif "exper" not in tok:
                raise ValueError("Unexpected token %s" % tok)
        if not k or not n:
            raise Exception("Invalid mashexp filename")
        try:
            k1 = basename(next(ifp).strip().split()[1])
        except StopIteration:
            print("StopIteration from file with path %s" % path)
            raise
        ret = {tuple(["_".join(sorted([k1, basename(i.split()[0])])), k, n]):
               float(i.strip().split()[1]) for i in ifp}
    return ret


def folder2dict(paths):
    if isinstance(paths, str):
        if os.path.isdir(paths):
            import glob
            paths = glob.glob(paths + "/*")
    ret = {}
    for path in paths:
        if "experim" not in path:
            continue
        ret.update(mash2dict(path))
    return ret


__all__ = ["mash2dict", "folder2dict", "get_flag"]
