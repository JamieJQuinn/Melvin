#!/usr/bin/env python3
"""Plotting functions for simulation variables"""
import sys
import numpy as np
from matplotlib import rc
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import argparse
import json

def main():
    """main function"""
    parser = argparse.ArgumentParser(description='Print variable contents')

    parser.add_argument('filename', help='file to open')
    parser.add_argument('--max_print_mode', type=int,
                        help='modes to print')
    parser.add_argument('--output',
                        help='Output file')
    parser.add_argument('--constants', help='constants file')
    # parser.add_argument('--plot_vorticity', action="store_true")
    args = parser.parse_args()

    constants_file = open(args.constants, "r")
    constants = json.load(constants_file)

    is_ddc = constants["isDoubleDiffusion"]

    n_modes = constants["nN"]
    n_z_gridpoints = constants["nZ"]
    aspect_ratio = constants["aspectRatio"]

    if "horizontalBoundaryConditions" not in constants:
        constants["horizontalBoundaryConditions"] = "impermeable"

    # Get data
    data = np.fromfile(args.filename, dtype=np.dtype(np.cdouble))
    temp = np.transpose(data[:n_modes*n_z_gridpoints]\
                        .reshape(n_modes, n_z_gridpoints))
    omg = np.transpose(data[1*n_modes*n_z_gridpoints:2*n_modes*n_z_gridpoints]\
                       .reshape(n_modes, n_z_gridpoints))
    psi = np.transpose(data[2*n_modes*n_z_gridpoints:3*n_modes*n_z_gridpoints]\
                       .reshape(n_modes, n_z_gridpoints))
    if is_ddc:
        varidx = 7
        xi = np.transpose(data[(varidx)*n_modes*n_z_gridpoints:(varidx+1)*n_modes*n_z_gridpoints]\
                           .reshape(n_modes, n_z_gridpoints))

    if constants["horizontalBoundaryConditions"] == "periodic":
        temp_actual = 2*(n_modes-1)*np.fft.irfft(temp)
        psi_actual = 2*(n_modes-1)*np.fft.irfft(psi)
        omg_actual = 2*(n_modes-1)*np.fft.irfft(omg)
    # else:
        # n_x_gridpoints = int(aspect_ratio*n_modes)
        # temp_actual = np.zeros((n_z_gridpoints, n_x_gridpoints))
        # x_axis = np.linspace(0, aspect_ratio, num=n_x_gridpoints)
        # cosine = np.cos(np.pi/aspect_ratio*np.outer(np.arange(n_modes), x_axis))
        # np.dot(np.real(temp), cosine, out=temp_actual)

    # Convert from spectral space to real space
    # psi_actual = np.zeros((n_z_gridpoints, n_x_gridpoints))
    # cosine = np.cos(np.pi/aspect_ratio*np.outer(np.arange(n_modes), x_axis))
    # sine = np.sin(np.pi/aspect_ratio*np.outer(np.arange(n_modes), x_axis))
    # np.dot(psi, sine, out=psi_actual)

    # if args.plot_vorticity:
        # omg_actual = np.zeros((n_z_gridpoints, n_x_gridpoints))
        # np.dot(omg, cosine, out=omg_actual)

    if is_ddc:
        xi_actual = np.zeros((n_z_gridpoints, n_x_gridpoints))
        np.dot(xi, cosine, out=xi_actual)

    N_PLOTS = 2

    fig, axes = plt.subplots(N_PLOTS, 3, figsize=(3*2*N_PLOTS*1.618, 2*N_PLOTS))

    # Plot
    z = np.linspace(0, 1, num=n_z_gridpoints)
    axes[0, 0].set_xticks([0, 1])
    axes[0, 0].plot(z, np.real(temp[:, 0]))
    axes[0, 0].set_xlabel(r"$z$")
    axes[0, 0].set_title(r"$T$")

    v_x = -(psi[2:, 0] - psi[:-2, 0])*(n_z_gridpoints-1)/2.0
    axes[0, 1].set_xticks([0, 1])
    axes[0, 1].plot(z[1:-1], np.real(v_x))
    axes[0, 1].set_xlabel(r"$z$")
    axes[0, 1].set_title(r"$\psi$")

    axes[0, 2].set_xticks([0, 1])
    axes[0, 2].plot(z, np.real(omg[:, 0]))
    axes[0, 2].set_xlabel(r"$z$")
    axes[0, 2].set_title(r"$\omega$")


    modes = range(0, n_modes)
    z_axis = np.linspace(0, 1, num=n_z_gridpoints)
    mode_grid, z_grid = np.meshgrid(modes, z_axis)
    # axes[1, 0].plot(modes, np.real(temp[n_z_gridpoints/2]))
    # axes[1, 1].plot(modes, np.real(psi[n_z_gridpoints/2]))
    # axes[1, 2].plot(modes, np.real(omg[n_z_gridpoints/2]))
    axes[1, 0].pcolormesh(mode_grid, z_grid, np.absolute(temp))
    axes[1, 1].pcolormesh(mode_grid, z_grid, np.absolute(psi))
    axes[1, 2].pcolormesh(mode_grid, z_grid, np.absolute(omg))
    axes[1, 0].set_xlim(0, 50)
    axes[1, 1].set_xlim(0, 50)
    axes[1, 2].set_xlim(0, 50)

    # x = np.linspace(0, constants["aspectRatio"], temp_actual.shape[1])
    # axes[2, 0].plot(x, temp_actual[n_z_gridpoints/2])
    # axes[2, 1].plot(x, psi_actual[n_z_gridpoints/2])
    # axes[2, 2].plot(x, omg_actual[n_z_gridpoints/2])

    if args.output:
        plt.savefig(args.output, bbox_inches='tight', pad_inches=0, dpi=200)
    else:
        plt.show()

main()
