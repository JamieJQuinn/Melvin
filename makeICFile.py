import numpy as np
import sys

nZ = 101
nN = 51

# Stored as temp|omg|psi contiguously
data = np.zeros(3*nN*nZ)
data[nZ*0+0:nZ*0+nZ] = np.linspace(1, 0, nZ)
for n in [1]:
	data[nZ*n+0:nZ*n+nZ] = 0.01*np.sin(np.pi*np.linspace(0, 1, nZ))

data.tofile(sys.argv[1])
