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
                        help='modes to initialise')
    parser.add_argument('--periodic', action='store_true',
                        help='enables periodic conditions')
    parser.add_argument('--salt_fingering', action='store_true',
                        help='enables salt fingering conditions')
    parser.add_argument('--combined_convection', action='store_true',
                        help='fully unstable convection in both gradients')
    parser.add_argument('--linear_stability', action='store_true',
                        help='sets up linear stability conditions')
    args = parser.parse_args()

    if args.salt_fingering or args.combined_convection:
        ddc = True

    n_modes = args.n_modes
    n_gridpoints = args.n_gridpoints
    if ddc:
        n_vars = 10
    else:
        n_vars = 7

    # Default to unstable temp gradient
    temp_grad = -1
    xi_grad = 1

    if args.salt_fingering:
        temp_grad = 1
        xi_grad = 1

    if args.combined_convection:
        temp_grad = -1
        xi_grad = 1

    # Stored as temp|omg|psi contiguously
    data = np.zeros(n_vars*n_modes*n_gridpoints)

    # Set up n=0 linear gradients
    if not args.periodic:
        mode = 0
        varidx = 0*n_modes*n_gridpoints # temperature
        data[varidx + n_gridpoints*(mode+0):varidx + n_gridpoints*(mode+1)] =\
                np.linspace(0, 1, n_gridpoints)[::temp_grad]

        if ddc:
            mode = 0
            varidx = 7*n_modes*n_gridpoints # xi
            data[varidx + n_gridpoints*(mode+0):varidx + n_gridpoints*(mode+1)] =\
                    np.linspace(0, 1, n_gridpoints)[::xi_grad]

    if args.linear_stability:
        # initialise all modes with a large amplitude
        modes = range(1, args.n_modes)
        amp = 1.0
    else:
        # initialise chosen modes with a small amplitude
        modes = args.modes
        amp = 0.01

    # Initialise modes
    for mode in modes:
        if not ddc:
            varidx = 0*n_modes*n_gridpoints # temperature
            data[varidx + n_gridpoints*(mode+0):varidx + n_gridpoints*(mode+1)] =\
                amp*np.sin(np.pi*np.linspace(0, 1, n_gridpoints))

        if ddc:
            varidx = 7*n_modes*n_gridpoints # xi
            data[varidx + n_gridpoints*(mode+0):varidx + n_gridpoints*(mode+1)] =\
                amp*np.sin(np.pi*np.linspace(0, 1, n_gridpoints))

    data.tofile(args.output)

main()
