#!/usr/bin/env python3
"""Prints mode, temperature, vorticity and stream function"""

import argparse
import numpy as np
import json

def main():
    """main function"""
    parser = argparse.ArgumentParser(description='Print variable contents')

    parser.add_argument('filenames', help='files to open', nargs='+')
    parser.add_argument('--max_print_mode', type=int,
                        help='modes to print')
    parser.add_argument('--constants', help='constants file')
    args = parser.parse_args()

    constants_file = open(args.constants, "r")
    constants = json.load(constants_file)

    is_ddc = constants["isDoubleDiffusion"]

    n_modes = constants["nN"]
    n_gridpoints = constants["nZ"]

    data = np.fromfile(args.filenames[0], dtype=np.dtype(np.double))
    temp = np.transpose(data[:n_modes*n_gridpoints]\
                        .reshape(n_modes, n_gridpoints))
    omega = np.transpose(data[n_modes*n_gridpoints:2*n_modes*n_gridpoints]\
                         .reshape(n_modes, n_gridpoints))
    psi = np.transpose(data[2*n_modes*n_gridpoints:3*n_modes*n_gridpoints]\
                       .reshape(n_modes, n_gridpoints))
    if(is_ddc):
        xi = np.transpose(data[7*n_modes*n_gridpoints:8*n_modes*n_gridpoints]\
                       .reshape(n_modes, n_gridpoints))

    for i in range(args.max_print_mode+1):
        if(not is_ddc):
            print(i, "%.5E" % temp[32, i], "%.5E" % omega[32, i], "%.5E" % psi[32, i])
        else:
            print(i, "%.5E" % temp[32, i], "%.5E" % omega[32, i], "%.5E" % psi[32, i], "%.5E" % xi[32, i])

main()
