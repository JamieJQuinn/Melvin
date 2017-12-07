#!/usr/bin/env python3
"""Prints mode, temperature, vorticity and stream function"""

import argparse
import numpy as np

def main():
    """main function"""
    parser = argparse.ArgumentParser(description='Print variable contents')

    parser.add_argument('filenames', help='files to open', nargs='+')
    parser.add_argument('--max_print_mode', type=int,
                        help='modes to print')
    parser.add_argument('--n_modes', type=int,
                        help='number of modes in simulation')
    parser.add_argument('--n_gridpoints', type=int,
                        help='number of gridpoints in simulation')
    args = parser.parse_args()

    n_modes = args.n_modes
    n_gridpoints = args.n_gridpoints
    data = np.fromfile(args.filenames[0], dtype=np.dtype(np.double))
    temp = np.transpose(data[:n_modes*n_gridpoints]\
                        .reshape(n_modes, n_gridpoints))
    omega = np.transpose(data[n_modes*n_gridpoints:2*n_modes*n_gridpoints]\
                         .reshape(n_modes, n_gridpoints))
    psi = np.transpose(data[2*n_modes*n_gridpoints:3*n_modes*n_gridpoints]\
                       .reshape(n_modes, n_gridpoints))

    for i in range(args.max_print_mode+1):
        print i, "%.5E" % temp[32, i], "%.5E" % omega[32, i], "%.5E" % psi[32, i]

main()
