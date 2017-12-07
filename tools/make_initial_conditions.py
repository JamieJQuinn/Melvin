#!/usr/bin/env python2
"""Create initial conditions file"""
import argparse
import numpy as np

def main():
    """main function"""
    parser = argparse.ArgumentParser(description='Create initial conditions file')

    parser.add_argument('--output', help='output file', required=True)
    parser.add_argument('--n_modes', type=int, required=True,
                        help='number of modes in simulation')
    parser.add_argument('--n_gridpoints', type=int, required=True,
                        help='number of gridpoints in simulation')
    parser.add_argument('--modes', type=int, nargs='+',
                        help='number of gridpoints in simulation')
    parser.add_argument('--periodic', action='store_true',
                        help='enables periodic conditions')
    args = parser.parse_args()

    n_modes = args.n_modes
    n_gridpoints = args.n_gridpoints

    # Stored as temp|omg|psi contiguously
    data = np.zeros(7*n_modes*n_gridpoints)
    if not args.periodic:
        mode = 0
        data[n_gridpoints*(mode+0):n_gridpoints*(mode+1)] = np.linspace(1, 0, n_gridpoints)

    for mode in args.modes:
        data[n_gridpoints*(mode+0):n_gridpoints*(mode+1)]\
                = 0.01*np.sin(np.pi*np.linspace(0, 1, n_gridpoints))

    data.tofile(args.output)

main()
