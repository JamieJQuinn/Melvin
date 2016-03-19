import numpy as np
import matplotlib.pyplot as plt
import sys

def main():
    for n, f in enumerate(sys.argv[2:]):
        ke = np.fromfile(f, dtype=np.dtype(np.double))
        plt.gca().set_yscale('log')
        plt.plot(ke, label=f.split('.')[0][17:])
    plt.show()

main()
