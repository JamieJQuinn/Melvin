#!/usr/bin/env python2
"""Plot kinetic energy"""
import argparse
import numpy as np
import matplotlib as mpl

def figsize(scale):
    """Used for latex calculations"""
    fig_width_pt = 512.14963
    inches_per_pt = 1.0/72.27
    golden_mean = (np.sqrt(5.0)-1.0)/2.0
    fig_width = fig_width_pt*inches_per_pt*scale
    fig_height = fig_width*golden_mean
    fig_size = [fig_width, fig_height]
    return fig_size

LATEX_CONFIG = {                  # setup matplotlib to use latex for output
    "pgf.texsystem": "pdflatex",    # change this if using xetex or lautex
    "text.usetex": True,            # use LaTeX to write all text
    "font.family": "serif",
    "font.serif": [],
    "font.sans-serif": [],
    "font.monospace": [],
    "axes.labelsize": 10,               # LaTeX default is 10pt font.
    "font.size": 10,
    "legend.fontsize": 8,       # Make the legend/label fonts a little smaller
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,

    "figure.figsize": figsize(1.0),     # default fig size of 0.9 textwidth
    "pgf.preamble": [
        r"\usepackage[utf8x]{inputenc}",    # use utf8 fonts becasue your computer can handle it :)
        r"\usepackage[T1]{fontenc}",        # plots will be generated using this preamble
        ]
    }
mpl.rcParams.update(LATEX_CONFIG)

import matplotlib.pyplot as plt

def main():
    """main function"""
    parser = argparse.ArgumentParser(description='Plot kinetic energy')

    parser.add_argument('filenames', nargs='+',
                        help='files to plot')
    parser.add_argument('--start_time', type=float, default=-1.,
                        help='average from here')
    parser.add_argument('--end_time', type=float, default=-1.,
                        help='average to here')
    parser.add_argument('--plot_against_mode', action='store_true',
                        help='plot average kinetic energy against mode')
    parser.add_argument('--output', help='output file')
    args = parser.parse_args()

    start_time = max(args.start_time, 0)
    end_time = args.end_time
    average = np.zeros(len(args.filenames))

    if not args.plot_against_mode:
        for filename in args.filenames:
            kinetic_energy = np.fromfile(filename, dtype=np.dtype(np.double))
            plt.gca().set_yscale('log')
            plt.plot(kinetic_energy, label="n="+filename.split('.')[0][17:])
        plt.xlabel(r"Time ($10^{-4}$)")
        plt.ylabel("Kinetic Energy")
        plt.tight_layout()
    else:
        for filename in args.filenames:
            mode = int(f.split('.')[0][17:])
            kinetic_energy = np.fromfile(f, dtype=np.dtype(np.double))
            if end_time <= start_time:
                kinetic_energy = kinetic_energy[start_time:]
            else:
                kinetic_energy = kinetic_energy[start_time:end_time]
            average[mode-1] = np.average(kinetic_energy)
        plt.gca().set_yscale('log')
        plt.plot(average)

    if args.output:
        plt.savefig(args.output)
    else:
        plt.show()

main()
