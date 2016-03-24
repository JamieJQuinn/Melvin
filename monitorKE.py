import numpy as np
import sys
import matplotlib
matplotlib.use('TKAgg')
import matplotlib.pyplot as plt
import time

def main():
	plt.ion()
	#plt.show()
	while True:
		for f in sys.argv[1:]:
			ke = np.fromfile(f, dtype=np.dtype(np.double))
			plt.gca().set_yscale('log')
			plt.plot(ke, label="n="+f.split('.')[0][17:])
		plt.legend(loc=0)
		plt.xlabel("Time (10^{-4})")
		plt.ylabel("Kinetic Energy")
		plt.tight_layout()
		plt.draw()
		time.sleep(5)

main()
