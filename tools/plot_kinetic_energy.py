#!/usr/bin/env python2
"""Plots kinetic energy"""
import sys
import numpy as np
import matplotlib.pyplot as plt

def main():
    """main function"""
    for filename in sys.argv[1:]:
        kinetic_energy = np.fromfile(filename, dtype=np.dtype(np.double))
        plt.gca().set_yscale('log')
        plt.plot(kinetic_energy, label="n="+filename.split('.')[0][17:])
    plt.xlabel("Time (10^{-4})")
    plt.ylabel("Kinetic Energy")
    plt.tight_layout()
    plt.show()

main()
