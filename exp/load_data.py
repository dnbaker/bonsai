#!/usr/bin/env python
import numpy as np


def load_data(path):
    return np.array([line.strip().split('\t')[2:] for
                     line in open(path) if line[0] != "#"],
                    dtype=np.double)


if __name__ == "__main__":
    import sys
    data = load_data(sys.argv[1])
    print("data: %s, %s" % (data, str(data.shape)))
