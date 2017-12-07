#!/usr/bin/env python2
"""Plotting functions for simulation variables"""
import sys
import numpy as np
import matplotlib.pyplot as plt


def main():
    """main function"""
    parser = argparse.ArgumentParser(description='Print variable contents')

    parser.add_argument('filenames', help='files to open', nargs='+')
    parser.add_argument('--max_print_mode', type=int,
                        help='modes to print')
    parser.add_argument('--n_modes', type=int,
                        help='number of modes in simulation')
    parser.add_argument('--n_z_gridpoints', type=int,
                        help='number of gridpoints in simulation')
    parser.add_argument('--aspect_ratio', type=int,
                        help='aspect ratio of sim')
    args = parser.parse_args()

    n_modes = args.n_modes
    n_z_gridpoints = args.n_gridpoints
    aspect_ratio = int(sys.argv[2])
    n_x_gridpoints = a*n_z_gridpoints
    x_axis = np.linspace(0, a, num=n_x_gridpoints)
    z_axis = np.linspace(0, 1, num=n_z_gridpoints)
    x_grid, z_grid = np.meshgrid(x_axis, z_axis)

    # Get data
    data = np.fromfile(sys.argv[1], dtype=np.dtype(np.double))
    temp = np.transpose(data[:n_modes*n_z_gridpoints]\
                        .reshape(n_modes, n_z_gridpoints))
    psi = np.transpose(data[2*n_modes*n_z_gridpoints:3*n_modes*n_z_gridpoints]\
                       .reshape(n_modes, n_z_gridpoints))
    #plt.gca().set_aspect(1.0)
    #plt.axis('off')
    #plt.gca().get_xaxis().set_visible(False)
    #plt.gca().get_yaxis().set_visible(False)
    #plt.plot(temp[:, 1])

    # Convert from spectral space to real space
    temp_actual = np.zeros((n_z_gridpoints, n_x_gridpoints))
    psi_actual = np.zeros((n_z_gridpoints, n_x_gridpoints))
    cosine = np.cos(np.pi/a*np.outer(np.arange(n_modes), x_axis))
    sine = np.sin(np.pi/a*np.outer(np.arange(n_modes), x_axis))
    np.dot(psi, cosine, out=psi_actual)
    np.dot(temp, sine, out=temp_actual)

    # Plot
    plt.subplot(1, 2, 1)
    plt.gca().set_aspect(1.0)
    plt.yticks([0, 1])
    plt.xticks(xrange(0, a+1))
    plt.tick_params(axis='x', bottom='off', top='off')
    plt.pcolormesh(x_grid, z_grid, temp_actual, cmap='coolwarm')
    plt.contour(x_grid, z_grid, temp_actual)
    #plt.savefig(sys.argv[1] +'tmp.png', bbox_inches='tight', pad_inches=0)
    plt.subplot(1, 2, 2)
    plt.gca().set_aspect(1.0)
    plt.gca().get_yaxis().set_visible(False)
    plt.xticks(xrange(0, a+1))
    plt.tick_params(axis='x', bottom='off', top='off')
    plt.gcf().tight_layout()
    plt.contour(x_grid, z_grid, psi_actual, 6, colors='k')
    plt.savefig(sys.argv[1] +'.png', bbox_inches='tight', pad_inches=0, dpi=200)
    #plt.show()

main()
