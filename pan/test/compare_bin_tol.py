#!/usr/bin/env python

"""
Compare two binary files and look for the maximum element-wise difference.
If the difference is greater than the tolerance print a message and exit.

"""
from __future__ import print_function
import os
import sys
import struct

import numpy as np


def read_file(filename):
    print("Reading: '%s' " % filename)
    with open(filename, "rb") as fh:
        numrows = struct.unpack('l', fh.read(8))[0]
        numcols = struct.unpack('l', fh.read(8))[0]
        data = np.zeros((numrows, numcols), dtype=np.float32)
        for j in range(numcols):
            for i in range(numrows):
                val = fh.read(4)
                data[i][j] = struct.unpack('f', val)[0]

    return data


def main(args):
    # read data
    data1 = read_file(args.file1)
    data2 = read_file(args.file2)
    
    # check sizes are the same
    if data1.shape[0] != data2.shape[0] or data1.shape[1] != data2.shape[1]:
        sys.stderr.write("Error: shapes of arrays do not match (%r != %r)\n" % (data1.shape, data2.shape))
        sys.exit(1)

    # tolerance
    tol = args.tol
    print("Tolerance is: %e" % tol)
    
    # compute max absolute difference between elements
    diff = np.subtract(data1, data2)
    absdiff = np.absolute(diff)
    maxdiff = absdiff.max()
    print("Max difference is: %e" % maxdiff)
    if maxdiff > args.tol:
        sys.stderr.write("Max error is greater than tolerance: %e > %e\n" % (maxdiff, tol))
        sys.exit(2)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Compare two binaries files using the given tolerance for the comparison.")
    parser.add_argument("file1", help="The first binary file")
    parser.add_argument("file2", help="The second binary file")
    parser.add_argument("-t", "--tol", type=float, default=1e-4, help="The tolerance to use for comparing elements (default=1e-4)")
    args = parser.parse_args()
    main(args)
