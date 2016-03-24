import numpy as np
import matplotlib.pyplot as plt
import sys

t0 = int(sys.argv[1])
t1 = int(sys.argv[2])
average = np.zeros(len(sys.argv[3:]))

for f in sys.argv[3:]:
	n = int(f.split('.')[0][17:])
	ke = np.fromfile(f, dtype=np.dtype(np.double))
	if t1==0:
		ke = ke[t0:]
	else:
		ke = ke[t0:t1]
	average[n-1] = np.average(ke)

plt.gca().set_yscale('log')
plt.plot(average)
plt.show()
