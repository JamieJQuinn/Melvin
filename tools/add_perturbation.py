#!/usr/bin/env python2
"""Adds perturbation to temperature variable"""

import argparse
import numpy as np

def main():
    """main function"""
    parser = argparse.ArgumentParser(description='Adds perturbation to temperature variable')

    parser.add_argument('--input', help='input file', required=True)
    parser.add_argument('--output', help='output file', required=True)
    parser.add_argument('--mode', type=int, required=True,
                        help='mode to add perturbation to')
    parser.add_argument('--n_modes', type=int, required=True,
                        help='number of modes in simulation')
    parser.add_argument('--n_gridpoints', type=int, required=True,
                        help='number of gridpoints in simulation')
    parser.add_argument('--amplitude', type=float, required=True,
                        help='amplitude of perturbation')
    args = parser.parse_args()

    n_modes = args.n_modes
    n_gridpoints = args.n_gridpoints
    mode = args.mode
    amp = args.amplitude
    data = np.fromfile(args.input, dtype=np.dtype(np.double))

    data[mode*n_gridpoints:(mode+1)*n_gridpoints]\
            += amp*np.sin(np.pi*np.linspace(0, 1, n_gridpoints))

    print "Adding " + str(amp) + " perturbation to mode "\
          + str(mode) + "; " + args.input + " -> " + args.output
    data.tofile(args.output)

main()
