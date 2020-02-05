#!/usr/bin/env python2
"""Create initial conditions file"""
import argparse
import numpy as np

def temperature_gradient(is_stable):
    if is_stable:
        return 1
    else:
        return -1

def xi_gradient(is_stable):
    if is_stable:
        return -1
    else:
        return 1

def set_temp(data, mode, values):
    data[0, mode, :] = values

def set_xi(data, mode, values):
    data[7, mode, :] = values

def set_temperature_background(data, background):
    set_temp(data, 0, background)

def set_xi_background(data, background):
    set_xi(data, 0, background)

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
    parser.add_argument('--fully_stable', action='store_true',
                        help='fully stable in both gradients')
    parser.add_argument('--thermal_rayleigh_taylor', action='store_true', default=False,
                        help='sets up thermal Rayleigh-Taylor conditions')
    parser.add_argument('--linear_stability', action='store_true',
                        help='sets up linear stability conditions')
    parser.add_argument('--step_profile', action='store_true',
                        help='uses step function instead of linear background')
    parser.add_argument('--amp', type=float,
                        help='amplitude of initial disturbance')

    args = parser.parse_args()

    ddc = args.salt_fingering or args.combined_convection or args.fully_stable

    n_modes = args.n_modes
    n_gridpoints = args.n_gridpoints
    if ddc:
        n_vars = 10
    else:
        n_vars = 7

    # Default to unstable temperature
    temp_grad = temperature_gradient(is_stable=False)
    xi_grad = xi_gradient(is_stable=True)

    if args.salt_fingering:
        temp_grad = temperature_gradient(is_stable=True)
        xi_grad = xi_gradient(is_stable=False)

    if args.combined_convection:
        temp_grad = temperature_gradient(is_stable=False)
        xi_grad = xi_gradient(is_stable=False)

    if args.fully_stable:
        temp_grad = temperature_gradient(is_stable=True)
        xi_grad = xi_gradient(is_stable=True)

    if args.linear_stability:
        # initialise all modes with a large amplitude
        modes = range(1, args.n_modes)
        amp = 1.0
    else:
        # initialise chosen modes with a small amplitude
        if args.modes:
            modes = args.modes
        else:
            modes = []
        if args.amp:
            amp = args.amp
        else:
            amp = 0.01

    # Stored as temp|omg|psi contiguously
    data = np.zeros((n_vars, n_modes, n_gridpoints))

    background = np.zeros(n_gridpoints)

    # Set up n=0 background
    if args.step_profile:
        # Step function
        background = np.zeros(n_gridpoints)
        background[int(n_gridpoints/2):] = 1 # Set upper half to 1
    elif not args.periodic:
        # Linear gradient
        background = np.linspace(0, 1, n_gridpoints)

    set_temperature_background(data, background[::temp_grad])
    if ddc:
        set_xi_background(data, background[::xi_grad])

    # Set up perturbations
    if args.step_profile:
        perturbation = np.zeros(n_gridpoints)
        perturbation[int(n_gridpoints/2)] = amp
    else:
        perturbation = amp*np.sin(np.pi*np.linspace(0, 1, n_gridpoints))

    for mode in modes:
        if not args.salt_fingering:
            set_temp(data, mode, perturbation)
        if ddc:
            set_xi(data, mode, perturbation)

    data.tofile(args.output)

main()
