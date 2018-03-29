#!/usr/bin/env python
import matplotlib.pyplot as plt
import sys
plt.switch_backend('agg')

if sys.version_info[0] < 2:
    raise Exception("Requires Python 3.")
SZ = 0.5

def main():
    data = [i.strip().split() for i in open(sys.argv[1])]
    of = sys.argv[2] if sys.argv[2:] else None
    if not of:
        raise Exception("Provide output suffix")
    n, k = list(map(int, data[0][6:8]))
    subsets = {key: [i for i in data if i[-1] == key] for key in set(i[-1] for i in data)}
    pkeys = list(set(tuple(sorted(i[:2])) for i in data))
    tmp = {}
    vd = {k: {} for k in subsets.keys()}
    for i in data:
        tmp[tuple(sorted(i[:2]))] = float(i[3])
        vd[i[-1]][tuple(sorted(i[:2]))] = float(i[2])
    vals = [tmp[key] for key in pkeys]
    orig = [vd['original'][key] for key in pkeys]
    imp = [vd['ertl_improved'][key] for key in pkeys]
    emle = [vd['ertl_mle'][key] for key in pkeys]
    ejmle = [vd['ertl_joint_mle'][key] for key in pkeys]
    fig, ax = plt.subplots()
    ax.scatter(vals, orig, c='b', s=0.5, marker='^', alpha=0.2)
    ax.scatter(vals, imp, c='g', s=0.5, marker='*', alpha=0.2)
    ax.scatter(vals, emle, c='y', s=0.5, marker='+', alpha=0.2)
    ax.scatter(vals, ejmle, c='r', s=0.5, marker='s', alpha=0.2)
    plt.savefig(of + ".png")
    

    #Path1  Path2   Approximate jaccard index   Exact jaccard index Absolute difference %difference from exact value    Sketch size Kmer size

if __name__ == "__main__":
    sys.exit(main())
