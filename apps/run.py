#!/usr/bin/env python3
"""Profile applications driver: `run.py setup`, then `run.py m1 m2 ...` or `run.py all`."""
import importlib
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bnsapps import common, mock  # noqa: E402

ANALYSES = ["m1", "m2", "m3", "m4", "m4cal", "m5", "m6", "r0"]


def setup():
    for name, gens in [("full", mock.genomes("full")), ("mock", mock.genomes("mock"))]:
        db, nodes, info = mock.build_db(name, gens)
        print(f"[setup] {name}: {db} {info}")


def main(argv):
    if not argv:
        print(__doc__)
        return 1
    todo = ANALYSES if argv == ["all"] else argv
    for a in todo:
        t0 = time.time()
        if a == "setup":
            setup()
        else:
            importlib.import_module(f"bnsapps.{a}").run()
        print(f"[run] {a}: {time.time() - t0:.0f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
