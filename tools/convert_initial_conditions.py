#!/usr/bin/env python2
"""Converts initial condition between resolutions"""
import argparse
import numpy as np

def main():
    """main function"""
    parser = argparse.ArgumentParser(description='Converts initial condition between resolutions')

    parser.add_argument('--input', help='input file', required=True)
    parser.add_argument('--output', help='output file', required=True)
    parser.add_argument('--n_modes_old', type=int, required=True,
                        help='input number of modes')
    parser.add_argument('--n_gridpoints_old', type=int, required=True,
                        help='input number of gridpoints')
    parser.add_argument('--n_modes_new', type=int, required=True,
                        help='output number of modes')
    parser.add_argument('--n_gridpoints_new', type=int, required=True,
                        help='output number of gridpoints')
    args = parser.parse_args()

    n_gridpoints_old = args.n_gridpoints_old
    n_modes_old = args.n_modes_old
    n_gridpoints_new = args.n_gridpoints_new
    n_modes_new = args.n_modes_new
    infile = args.input
    outfile = args.output

    if (n_gridpoints_old > n_gridpoints_new) or (n_modes_old > n_modes_new):
        print "Can't convert to a lower resolution."
        exit(-1)

    print "Converting " + infile + " (" + str(n_gridpoints_old) + ", "\
            + str(n_modes_old) + ") to "\
            + outfile + " (" + str(n_gridpoints_new) + ", " + str(n_modes_new) + ")"

    z_ratio = float(n_gridpoints_new - 1)/(n_gridpoints_old - 1)

    data_old = np.fromfile(infile, dtype=np.dtype(np.double))
    data_new = np.zeros(7*n_gridpoints_new*n_modes_new)

    for i in xrange(7):
        for mode in xrange(n_modes_old):
            # interpolate
            data_new[(i*n_modes_new+mod)*n_gridpoints_new:\
                     (i*n_modes_new+mode+1)*n_gridpoints_new] = \
                np.interp(\
                np.linspace(0, 1, num=n_gridpoints_new), \
                np.linspace(0, 1, num=n_gridpoints_old), \
                data_old[(i*n_modes_old+mode)*n_gridpoints_old:\
                         (i*n_modes_old+mode+1)*n_gridpoints_old]\
                )


    data_new.tofile(outfile)

main()
