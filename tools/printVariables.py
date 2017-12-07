#!/usr/bin/env python3

import argparse
import numpy as np

def main():
    parser = argparse.ArgumentParser(description='Print variable contents')

    parser.add_argument('filenames', help='files to open', nargs='+')
    parser.add_argument('--max_print_mode', type=int,
                        help='modes to print')
    parser.add_argument('--n_modes', type=int,
                        help='number of modes in simulation')
    parser.add_argument('--n_gridpoints', type=int,
                        help='number of gridpoints in simulation')
    args = parser.parse_args()

    Nn = int(args.n_modes)
    Nz = int(args.n_gridpoints)
    data = np.fromfile(args.filenames[0], dtype=np.dtype(np.double))
    temp = np.transpose(data[:Nn*Nz].reshape(Nn, Nz))
    omega = np.transpose(data[Nn*Nz:2*Nn*Nz].reshape(Nn, Nz))
    psi = np.transpose(data[2*Nn*Nz:3*Nn*Nz].reshape(Nn, Nz))

    for i in range(args.max_print_mode+1):
        print i, "%.5E" % temp[32, i], "%.5E" % omega[32, i], "%.5E" % psi[32, i]

main()
