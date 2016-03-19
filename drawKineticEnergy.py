import numpy as np
import matplotlib as mpl
import sys
def figsize(scale):
    fig_width_pt = 512.14963
    inches_per_pt = 1.0/72.27
    golden_mean = (np.sqrt(5.0)-1.0)/2.0            
    fig_width = fig_width_pt*inches_per_pt*scale    
    fig_height = fig_width*golden_mean             
    fig_size = [fig_width,fig_height]
    return fig_size

pgf_with_latex = {                      # setup matplotlib to use latex for output
    "pgf.texsystem": "pdflatex",        # change this if using xetex or lautex
    "text.usetex": True,                # use LaTeX to write all text
    "font.family": "serif",
    "font.serif": [],                   # blank entries should cause plots to inherit fonts from the document
    "font.sans-serif": [],
    "font.monospace": [],
    "axes.labelsize": 10,               # LaTeX default is 10pt font.
    "text.fontsize": 10,
    "legend.fontsize": 8,               # Make the legend/label fonts a little smaller
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "figure.figsize": figsize(0.6),     # default fig size of 0.9 textwidth
    "pgf.preamble": [
        r"\usepackage[utf8x]{inputenc}",    # use utf8 fonts becasue your computer can handle it :)
        r"\usepackage[T1]{fontenc}",        # plots will be generated using this preamble
        ]
    }
mpl.rcParams.update(pgf_with_latex)

import matplotlib.pyplot as plt

def main():
    for f in sys.argv[1:]:
        ke = np.fromfile(f, dtype=np.dtype(np.double))
        plt.gca().set_yscale('log')
        plt.plot(ke, label="n="+f.split('.')[0][17:])
    plt.legend(loc=0)
    plt.xlabel(r"Time ($10^{-4}$)")
    plt.ylabel("Kinetic Energy")
    plt.xlim(0, 2000)
    plt.ylim(1e-3, 1e7)
    plt.tight_layout()
    plt.savefig('{}.pdf'.format("KineticEnergy"))
    plt.show()

main()
