import numpy as np
import matplotlib.pyplot as plt
import sys

t0 = int(sys.argv[1])
average = np.zeros(len(sys.argv[2:]))

for n, f in enumerate(sys.argv[2:]):
	ke = np.fromfile(f, dtype=np.dtype(np.double))
	ke = ke[t0:]
	average[n] = np.average(ke)

plt.gca().set_yscale('log')
plt.plot(average)
plt.show()
