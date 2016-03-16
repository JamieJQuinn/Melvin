import numpy as np
import matplotlib.pyplot as plt
import sys

def main():
    for f in sys.argv[1:]:
        ke = np.fromfile(f, dtype=np.dtype(np.double))
        plt.gca().set_yscale('log')
        plt.plot(ke)
    plt.show()

main()
