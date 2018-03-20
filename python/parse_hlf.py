import sys
import numpy as np


def main():
    with open(sys.argv[1]) if sys.argv[1:] else sys.stdin as f:
        data = list(f)
        headers = data[0][1:].split('\t')
        lines = np.array([list(map(float, line.split('\t'))) for
                          line in data[1:]], dtype=np.double)
        nums = set(lines[:,-3].astype(np.int))
        for num in nums:
            subset = lines[lines[:,-3] == num]
            meanerr = np.mean(subset[:,-5])
            print("meanerr: %lf. num: %i" % (meanerr, num))


if __name__ == "__main__":
    sys.exit(main())
