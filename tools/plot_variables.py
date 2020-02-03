#!/usr/bin/env python3
"""Plotting functions for simulation variables"""
import sys
import numpy as np
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
    args = parser.parse_args()

    constants_file = open(args.constants, "r")
    constants = json.load(constants_file)

    is_ddc = constants["isDoubleDiffusion"]

    n_modes = constants["nN"]
    n_z_gridpoints = constants["nZ"]
    aspect_ratio = constants["aspectRatio"]

    n_x_gridpoints = int(aspect_ratio*n_z_gridpoints)
    x_axis = np.linspace(0, aspect_ratio, num=n_x_gridpoints)
    z_axis = np.linspace(0, 1, num=n_z_gridpoints)
    x_grid, z_grid = np.meshgrid(x_axis, z_axis)

    # Get data
    data = np.fromfile(args.filename, dtype=np.dtype(np.float64))
    temp = np.transpose(data[:n_modes*n_z_gridpoints]\
                        .reshape(n_modes, n_z_gridpoints))
    psi = np.transpose(data[2*n_modes*n_z_gridpoints:3*n_modes*n_z_gridpoints]\
                       .reshape(n_modes, n_z_gridpoints))
    if is_ddc:
        varidx = 7
        xi = np.transpose(data[(varidx)*n_modes*n_z_gridpoints:(varidx+1)*n_modes*n_z_gridpoints]\
                           .reshape(n_modes, n_z_gridpoints))
    #plt.gca().set_aspect(1.0)
    #plt.axis('off')
    #plt.gca().get_xaxis().set_visible(False)
    #plt.gca().get_yaxis().set_visible(False)
    #plt.plot(temp[:, 1])

    # Convert from spectral space to real space
    temp_actual = np.zeros((n_z_gridpoints, n_x_gridpoints))
    psi_actual = np.zeros((n_z_gridpoints, n_x_gridpoints))
    cosine = np.cos(np.pi/aspect_ratio*np.outer(np.arange(n_modes), x_axis))
    sine = np.sin(np.pi/aspect_ratio*np.outer(np.arange(n_modes), x_axis))
    np.dot(psi, sine, out=psi_actual)
    np.dot(temp, cosine, out=temp_actual)

    if is_ddc:
        xi_actual = np.zeros((n_z_gridpoints, n_x_gridpoints))
        np.dot(xi, cosine, out=xi_actual)

    # Plot
    plt.subplot(1, 3, 1)
    plt.gca().set_aspect(1.0)
    plt.yticks([0, 1])
    plt.xticks(range(0, int(aspect_ratio+1)))
    plt.tick_params(axis='x', bottom=False, top=False)
    plt.pcolormesh(x_grid, z_grid, temp_actual, cmap='coolwarm')
    # plt.contour(x_grid, z_grid, temp_actual)
    #plt.savefig(sys.argv[1] +'tmp.png', bbox_inches='tight', pad_inches=0)

    if np.any(psi_actual):
        plt.subplot(1, 3, 2)
        plt.gca().set_aspect(1.0)
        plt.gca().get_yaxis().set_visible(False)
        plt.xticks(range(0, int(aspect_ratio+1)))
        plt.tick_params(axis='x', bottom=False, top=False)
        plt.gcf().tight_layout()
        plt.contour(x_grid, z_grid, psi_actual, 6, colors='k')

    if is_ddc:
        plt.subplot(1, 3, 3)
        plt.gca().set_aspect(1.0)
        plt.yticks([0, 1])
        plt.xticks(range(0, int(aspect_ratio+1)))
        plt.tick_params(axis='x', bottom=False, top=False)
        plt.pcolormesh(x_grid, z_grid, xi_actual, cmap='Blues')

    if args.output:
        plt.savefig(args.output, bbox_inches='tight', pad_inches=0, dpi=200)
    else:
        plt.show()

main()
