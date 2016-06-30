#!/usr/bin/env python

"""
Calculate differences in peaks between two cellsim output files.

"""
from __future__ import print_function
import sys
import math

from scipy import signal

import cs


def check_peak_times(curve1, maxima1, curve2, maxima2):
    """Compare the times of the maxima."""
    num_max1 = len(maxima1)
    num_max2 = len(maxima2)
    num_max = min(num_max1, num_max2)
    for i in range(num_max):
        step1 = maxima1[i]
        step2 = maxima2[i]
        # allow one time step offset (better as percentage?)
        if math.fabs(step1 - step2) > 1:
            sys.stderr.write("Error: peak {0} location differs ({1} - {2})\n".format(i, step1, step2))
            sys.exit(2)


def check_peak_values(curve1, maxima1, curve2, maxima2, args):
    """Compare the values of the maxima."""
    num_max1 = len(maxima1)
    num_max2 = len(maxima2)
    num_max = min(num_max1, num_max2)
    for i in range(num_max):
        val1 = curve1[maxima1[i]]
        val2 = curve2[maxima2[i]]
        percdiff = math.fabs((val1 - val2) / val2) * 100.0
        if percdiff > args.peaktol:
            sys.stderr.write("Error: peak {0} value differs ({1} vs {2} => {3} %)\n".format(i, val1, val2, percdiff))
            sys.exit(3)


def main(args):
    # read data
    print("Reading new result file: {0}".format(args.filenew))
    data1 = cs.get_data(args.filenew)
    print("Reading reference result file: {0}".format(args.fileref))
    data2 = cs.get_data(args.fileref)
    
    # check sizes are the same
    if data1.shape[0] != data2.shape[0] or data1.shape[1] != data2.shape[1]:
        sys.stderr.write("Error: shapes of arrays do not match (%r != %r)\n" % (data1.shape, data2.shape))
        sys.exit(1)
    
    # compare each curve
    for i in range(data1.shape[0]):
        print("Analysing curve {0}".format(i))
        
        # get the curves
        curve1 = data1[i][:]
        curve2 = data2[i][:]
        
        # finding maxima
        maxima1 = signal.argrelmax(curve1)[0]
        print("  Curve 1 maxima (%d) at steps: %r" % (len(maxima1), list(maxima1)))
        maxima2 = signal.argrelmax(curve2)[0]
        print("  Curve 2 maxima (%d) at steps: %r" % (len(maxima2), list(maxima2)))
        
        # check maxima occur at the same place
        print("  Analysing peak times...")
        check_peak_times(curve1, maxima1, curve2, maxima2)
        print("    Peak times OK")
        
        # check maxima values match
        print("  Analysing peak values (tolerance is {0} %)...".format(args.peaktol))
        check_peak_values(curve1, maxima1, curve2, maxima2, args)
        print("    Peak values OK")
    
    print("Peaks match!")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Compare peak locations and values in two cellsim files.")
    parser.add_argument("filenew", help="The new binary file")
    parser.add_argument("fileref", help="The reference binary file")
    parser.add_argument("-p", "--peaktol", type=float, default=1e-2, help="The percentage tolerance used for comparing maxima values (default=0.01%%)")
    args = parser.parse_args()
    main(args)
