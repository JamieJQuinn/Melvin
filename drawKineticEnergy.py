import numpy as np
import matplotlib.pyplot as plt
import sys

def main():
    ke = np.fromfile(sys.argv[1], dtype=np.dtype(np.double))
    plt.gca().set_yscale('log')
    plt.plot(ke)
    plt.show()

main()
