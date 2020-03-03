#!/usr/bin/env python
from operator import or_
from functools import reduce


def gen_val(i, nnucs):
    a, c, g, t = [0] * 4
    while nnucs:
        val = i & 0x3
        if val == 0:
            a += 1
        elif val == 1:
            c += 1
        elif val == 2:
            g += 1
        else:
            t += 1
        i >>= 2
        nnucs -= 1
    return reduce(or_,
                  ((i << (x << 3)) for
                   i, x in zip((a, c, g, t),
                               range(3, -1, -1))))


def printval(i, nnucs):
    a, c, g, t = i >> 24 & 0xFF, i >> 16 & 0xFF, i >> 8 & 0xFF, i & 0xFF
    print("value: %i. nnucs: %i. a: %i, c: %i, g: %i, t: %i." %
          (i, nnucs, a, c, g, t))


def generate_table(nnucs):
    vals = [0] * (1 << (nnucs << 1))
    for i in range(len(vals)):
        vals[i] = gen_val(i, nnucs)
    return vals


def table2str(vals):
    sublists = [vals[i << 4:(i+1) << 4] for i in range(len(vals) >> 4)]
    nnucs = 0
    tmp = len(vals) >> 2
    while tmp:
        nnucs += 1
        tmp >>= 2
    if not sublists:
        sublists = [vals]
    substrs = ["    %s," % ", ".join(map(str, sublist)) for
               sublist in sublists]
    return "static const u32 lut%i [] {\n%s\n};" % (
        nnucs, "\n".join(substrs) + "\n};")


if __name__ == "__main__":
    from sys import argv
    nnucs = int(argv[1]) if argv[1:] else 4
    vals = generate_table(nnucs)
    print(table2str(vals))
    if argv[2:]:
        {printval(val, nnucs) for val in vals}
