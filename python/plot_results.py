#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import sys
plt.switch_backend('agg')

if sys.version_info[0] < 2:
    raise Exception("Requires Python 3.")


def main():
    data = [i.strip().split() for i in open(sys.argv[1])]
    of = sys.argv[2] if sys.argv[2:] else None
    if not of:
        raise Exception("Provide output suffix")
    n, k = list(map(int, data[0][6:8]))
    subsets = {key: [i for i in data if i[-1] == key] for
               key in set(i[-1] for i in data)}
    pkeys = list(set(tuple(sorted(i[:2])) for i in data))
    tmp = {}
    vd = {k: {} for k in subsets.keys()}
    for i in data:
        tmp[tuple(sorted(i[:2]))] = float(i[3])
        vd[i[-1]][tuple(sorted(i[:2]))] = float(i[2])
    vals = np.array([tmp[key] for key in pkeys])
    orig = np.array([vd['original'][key] for key in pkeys])
    imp = np.array([vd['ertl_improved'][key] for key in pkeys])
    emle = np.array([vd['ertl_mle'][key] for key in pkeys])
    ejmle = np.array([vd['ertl_joint_mle'][key] for key in pkeys])
    fig, ax = plt.subplots()
    ax.scatter(vals, orig - vals, c='b', s=1.0, marker='^', alpha=0.4)
    ax.scatter(vals, imp - vals, c='g', s=1.0, marker='*', alpha=0.4)
    ax.scatter(vals, emle - vals, c='y', s=1.0, marker='+', alpha=0.4)
    ax.scatter(vals, ejmle - vals, c='r', s=1.0, marker='s', alpha=0.4)
    plt.savefig(of + ".png", dpi=1000)
    mask = vals > 0.01
    morig = orig[mask]
    mvals = vals[mask]
    memle = emle[mask]
    mimp = imp[mask]
    mejmle = ejmle[mask]
    ax.clear()
    ax.scatter(mvals, orig[mask] - mvals, c='b', s=1.0, marker='^', alpha=0.4)
    ax.scatter(mvals, imp[mask] - mvals, c='g', s=1.0, marker='*', alpha=0.4)
    ax.scatter(mvals, emle[mask] - mvals, c='y', s=1.0, marker='+', alpha=0.4)
    ax.scatter(mvals, ejmle[mask] - mvals, c='r', s=1.0, marker='s', alpha=0.4)
    plt.savefig(of + "gt.01.png", dpi=1000)
    ax.clear()
    ax.scatter(orig - vals, vals, c='b', s=1.0, marker='^', alpha=0.4)
    ax.scatter(imp - vals, vals, c='g', s=1.0, marker='*', alpha=0.4)
    ax.scatter(emle - vals, vals, c='y', s=1.0, marker='+', alpha=0.4)
    ax.scatter(ejmle - vals, vals, c='r', s=1.0, marker='s', alpha=0.4)
    plt.savefig(of + ".rev.png", dpi=1000)
    ax.scatter(morig - mvals, mvals, c='b', s=1.0, marker='^', alpha=0.4)
    ax.scatter(mimp - mvals, mvals, c='g', s=1.0, marker='*', alpha=0.4)
    ax.scatter(memle - mvals, mvals, c='y', s=1.0, marker='+', alpha=0.4)
    ax.scatter(mejmle - mvals, mvals, c='r', s=1.0, marker='s', alpha=0.4)
    plt.savefig(of + ".rev.gt.01.png", dpi=1000)
    sses = {"orig": np.dot(orig - vals, orig - vals),
            "imp": np.dot(imp - vals, imp - vals),
            "emle": np.dot(emle - vals, emle - vals),
            "ejmle": np.dot(ejmle - vals, ejmle - vals)}
    print(f"sum of squared errors: {sses}")
    msses = {"orig": np.dot(morig - mvals, morig - mvals),
             "imp": np.dot(mimp - mvals, mimp - mvals),
             "emle": np.dot(memle - mvals, memle - mvals),
             "ejmle": np.dot(mejmle - mvals, mejmle - mvals)}
    print(f"sum of squared errors, masked > 0.01: {msses}")
    biases = {"orig": np.sum(orig - vals),
              "imp": np.sum(imp - vals),
              "emle": np.sum(emle - vals),
              "ejmle": np.sum(ejmle - vals)}
    print(f"biases: {biases}")
    mbiases = {"orig": np.sum(morig - mvals),
               "imp": np.sum(mimp - mvals),
               "emle": np.sum(memle - mvals),
               "ejmle": np.sum(mejmle - mvals)}
    print(f"biases, masked > 0.01: {mbiases}")


if __name__ == "__main__":
    sys.exit(main())
