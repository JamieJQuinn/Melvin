#!/usr/bin/env python3
"""Plotting functions for simulation variables"""
import numpy as np
from matplotlib import rc
rc('text', usetex=True)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import argparse
import json

from helper_functions import get_spatial_data

def set_plot_defaults(axis, aspect_ratio):
    axis.set_aspect(1.0)
    axis.set_yticks([0.0, 1.0])
    axis.set_xticks(range(0, int(aspect_ratio+1)))
    axis.set_ylabel(r"$z$")
    # axis.set_xlabel(r"$x$")

def main():
    """main function"""
    parser = argparse.ArgumentParser(description='Print variable contents')

    parser.add_argument('filename', help='file to open')
    parser.add_argument('--output',
                        help='Output file')
    parser.add_argument('--constants', help='constants file')
    parser.add_argument('--plot_vorticity', action="store_true")
    parser.add_argument('--plot_streamfunction', action="store_true")
    parser.add_argument('--colourbars', action="store_true")
    args = parser.parse_args()

    constants_file = open(args.constants, "r")
    constants = json.load(constants_file)

    is_ddc = constants["isDoubleDiffusion"]
    n_modes = constants["nN"]
    n_z_gridpoints = constants["nZ"]
    aspect_ratio = constants["aspectRatio"]

    if "horizontalBoundaryConditions" not in constants:
        constants["horizontalBoundaryConditions"] = "impermeable"

    n_plots = 1
    if args.plot_streamfunction:
        n_plots += 1
    if args.plot_vorticity:
        n_plots += 1
    if is_ddc:
        n_plots += 1

    # Get data
    data = np.fromfile(args.filename, dtype=np.dtype(np.cdouble))
    temp = get_spatial_data(data, 0, constants, False)
    if args.plot_streamfunction:
        psi = get_spatial_data(data, 2, constants, True)
    if args.plot_vorticity:
        omg = get_spatial_data(data, 1, constants, True)
    if is_ddc:
        xi = get_spatial_data(data, 7, constants, False)

    # Plot
    n_x_gridpoints = temp.shape[1]
    x_axis = np.linspace(0, aspect_ratio, num=n_x_gridpoints)
    z_axis = np.linspace(0, 1, num=n_z_gridpoints)
    x_grid, z_grid = np.meshgrid(x_axis, z_axis)

    _fig, axes = plt.subplots(n_plots, sharex=True)

    plot_index = 0
    set_plot_defaults(axes[plot_index], aspect_ratio)
    axes[plot_index].pcolormesh(x_grid, z_grid, temp, shading='auto', cmap='RdBu_r')
    axes[plot_index].set_title(r'$T$')

    if args.plot_streamfunction:
        plot_index += 1
        set_plot_defaults(axes[plot_index], aspect_ratio)
        axes[plot_index].contour(x_grid, z_grid, psi, 20, colors='k', linewidths=0.25)
        axes[plot_index].set_title(r'$\psi$')

    if is_ddc:
        plot_index += 1
        set_plot_defaults(axes[plot_index], aspect_ratio)
        axes[plot_index].set_title(r"$\xi$")
        axes[plot_index].pcolormesh(x_grid, z_grid, xi, shading='auto', cmap='Blues')

    if args.plot_vorticity:
        plot_index += 1
        set_plot_defaults(axes[plot_index], aspect_ratio)
        axes[plot_index].set_title(r"$\omega$")
        axes[plot_index].pcolormesh(x_grid, z_grid, omg, shading='auto', cmap='RdGy_r')

    axes[plot_index].set_xlabel(r"$x$")

    # Save
    if args.output:
        plt.savefig(args.output, bbox_inches='tight', pad_inches=0, dpi=200)
    else:
        plt.show()

main()
