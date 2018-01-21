#!/usr/bin/env python

if __name__ == "__main__":
    from sys import stderr
    w = stderr.write
    w("static const unsigned uint8lut [] {\n")
    w("    " + ", ".join(map(str, (bin(x).count("1") for x in range(256)))))
    w("\n};\n")
