#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import sys
plt.switch_backend('agg')

if sys.version_info[0] < 2:
    raise Exception("Requires Python 3.")
if sys.version_info[1] < 6:
    raise Exception("Requires Python >= 3.6")


def get_args():
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("-o", "--outsuffix", default="ashcmp")
    p.add_argument("folder", help="data folder path")
    p.add_argument("-p", "--threads", type=int, default=8)
    return p.parse_args()


def folder2paths(path):
    import glob
    import os
    fdpaths = glob.glob(path + "/flashdans*tsv")
    if not fdpaths:
        raise Exception("The files are NOT in the computer.")
    combs = set(tuple(map(int, i.split(".")[1:3])) for i in fdpaths)
    mpaths = glob.glob(path + "/mash.[0-9]*tsv")
    if not mpaths:
        raise Exception("The files are NOT in the computer.")
    mcombs = set(tuple(map(int, i.split(".")[1:3])) for i in mpaths)
    retcombs = combs & mcombs
    ret = list((f"{path}/flashdans.{c[0]}.{c[1]}.tsv",
                f"{path}/mash.{c[0]}.{c[1]}.tsv") for c in retcombs)
    for p1, p2 in ret:
        assert os.path.isfile(p1), f"p1 is not a file {p1}"
        assert os.path.isfile(p2), f"p2 is not a file {p2}"
    return ret


def main():
    import multiprocessing as mp
    args = get_args()
    Spooool = mp.Pool(args.threads)
    with open("summary." + args.outsuffix + ".txt", "w") as sfp:
        argsets = ((f, m, args.outsuffix) for f, m in folder2paths(args.folder))
        for outstr in Spooool.map(submit_subplot, argsets):
            sfp.write(outstr)


def get_nk(path):
    toks = path.split(".")
    ret = (int(toks[1]), int(toks[2]))
    return ret


def submit_subplot(tup):
    f, m, suffix = tup
    subplot(f, m, suffix)


def subplot(flashdata, mashdata, outsuffix):
    k, n = get_nk(flashdata)
    of = f"n{n}.k{k}.{outsuffix}"
    data = [i.strip().split() for i in open(flashdata)]
    mdata = [i.strip().split() for i in open(mashdata)]
    k_, n_ = list(map(int, data[0][6:8]))
    assert n == n_
    assert k == k_
    nkstr = f"[n:{n}|k:{k}]"
    subsets = {key: [i for i in data if i[-1] == key] for
               key in set(i[-1] for i in data)}
    pkeys = list(set(tuple(sorted(i[:2])) for i in data))
    tmp = {}
    vd = {k: {} for k in subsets.keys()}
    vd['mash'] = {}
    for i in data:
        tmp[tuple(sorted(i[:2]))] = float(i[3])
        vd[i[-1]][tuple(sorted(i[:2]))] = float(i[2])
    for i in mdata:
        vd['mash'][tuple(sorted(i[:2]))] = float(i[5])
    vals = np.array([tmp[key] for key in pkeys])
    orig = np.array([vd['original'][key] for key in pkeys])
    imp = np.array([vd['ertl_improved'][key] for key in pkeys])
    emle = np.array([vd['ertl_mle'][key] for key in pkeys])
    ejmle = np.array([vd['ertl_joint_mle'][key] for key in pkeys])
    mash = np.array([vd['mash'][key] for key in pkeys])
    fig, ax = plt.subplots()
    ax.scatter(vals, orig - vals, c='b', s=1.0, marker='^', alpha=0.4)
    ax.scatter(vals, imp - vals, c='g', s=1.0, marker='*', alpha=0.4)
    ax.scatter(vals, emle - vals, c='y', s=1.0, marker='+', alpha=0.4)
    ax.scatter(vals, ejmle - vals, c='r', s=1.0, marker='s', alpha=0.4)
    ax.scatter(vals, mash - vals, c='peru', s=1.0, marker='8', alpha=0.4)
    plt.savefig(of + ".png", dpi=1000)
    mask = vals > 0.01
    morig = orig[mask]
    mvals = vals[mask]
    memle = emle[mask]
    mimp = imp[mask]
    mejmle = ejmle[mask]
    mmash = mash[mask]
    ax.clear()
    x.scatter(mvals, morig - mvals, c='b', s=1.0, marker='^', alpha=0.4)
    x.scatter(mvals, mimp - mvals, c='g', s=1.0, marker='*', alpha=0.4)
    x.scatter(mvals, memle - mvals, c='y', s=1.0, marker='+', alpha=0.4)
    ax.scatter(mvals, mejmle - mvals, c='r', s=1.0, marker='s', alpha=0.4)
    ax.scatter(mvals, mmash - mvals, c='peru', s=1.0, marker='8', alpha=0.4)
    plt.savefig(of + "gt.01.png", dpi=1000)
    ax.clear()
    ax.scatter(orig - vals, vals, c='b', s=1.0, marker='^', alpha=0.4)
    ax.scatter(imp - vals, vals, c='g', s=1.0, marker='*', alpha=0.4)
    ax.scatter(emle - vals, vals, c='y', s=1.0, marker='+', alpha=0.4)
    ax.scatter(ejmle - vals, vals, c='r', s=1.0, marker='s', alpha=0.4)
    ax.scatter(mash - vals, vals, c='peru', s=1.0, marker='8', alpha=0.4)
    plt.savefig(of + ".rev.png", dpi=1000)
    ax.scatter(morig - mvals, mvals, c='b', s=1.0, marker='^', alpha=0.4)
    ax.scatter(mimp - mvals, mvals, c='g', s=1.0, marker='*', alpha=0.4)
    ax.scatter(memle - mvals, mvals, c='y', s=1.0, marker='+', alpha=0.4)
    ax.scatter(mejmle - mvals, mvals, c='r', s=1.0, marker='s', alpha=0.4)
    ax.scatter(mmash - mvals, mvals, c='pertu', s=1.0, marker='8', alpha=0.4)
    plt.savefig(of + ".rev.gt.01.png", dpi=1000)
    sses = {"orig": np.dot(orig - vals, orig - vals),
            "imp": np.dot(imp - vals, imp - vals),
            "emle": np.dot(emle - vals, emle - vals),
            "ejmle": np.dot(ejmle - vals, ejmle - vals),
            "mash": np.dot(mash - vals, mash - vals)}
    ret = f"{nkstr} sum of squared errors: {sses}\n"
    msses = {"orig": np.dot(morig - mvals, morig - mvals),
             "imp": np.dot(mimp - mvals, mimp - mvals),
             "emle": np.dot(memle - mvals, memle - mvals),
             "ejmle": np.dot(mejmle - mvals, mejmle - mvals),
             "mash": np.dot(mmash - mvals, mmash - mvals)}
    ret += f"{nkstr} sum of squared errors, masked > 0.01: {msses}\n"
    biases = {"orig": np.sum(orig - vals),
              "imp": np.sum(imp - vals),
              "emle": np.sum(emle - vals),
              "ejmle": np.sum(ejmle - vals),
              "mash": np.sum(ejmle - vals)}
    ret += f"{nkstr} biases: {biases}"
    mbiases = {"orig": np.sum(morig - mvals),
               "imp": np.sum(mimp - mvals),
               "emle": np.sum(memle - mvals),
               "ejmle": np.sum(mejmle - mvals),
               "mash": np.sum(mash - mvals)}
    ret += f"{nkstr} biases, masked > 0.01: {mbiases}"
    errors = {"orig": np.sum(np.abs(orig - vals)),
              "imp": np.sum(np.abs(imp - vals)),
              "emle": np.sum(np.abs(emle - vals)),
              "ejmle": np.sum(np.abs(ejmle - vals)),
              "mash": np.sum(np.abs(ejmle - vals))}
    ret += f"{nkstr} errors: {errors}"
    merrors = {"orig": np.sum(np.abs(morig - mvals)),
               "imp": np.sum(np.abs(mimp - mvals)),
               "emle": np.sum(np.abs(memle - mvals)),
               "ejmle": np.sum(np.abs(mejmle - mvals)),
               "mash": np.sum(np.abs(mash - mvals))}
    ret += f"{nkstr} errors, masked > 0.01: {merrors}"


if __name__ == "__main__":
    sys.exit(main())
