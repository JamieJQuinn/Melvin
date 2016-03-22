import numpy as np
import sys

nZ = int(sys.argv[1])
nN = int(sys.argv[2])
icModes = [int(arg) for arg in sys.argv[3:]]

# Stored as temp|omg|psi contiguously
data = np.zeros(7*nN*nZ)
data[nZ*0+0:nZ*0+nZ] = np.linspace(1, 0, nZ)
for n in icModes:
	data[nZ*n+0:nZ*n+nZ] = 0.01*np.sin(np.pi*np.linspace(0, 1, nZ))

data.tofile("ICn"+'_'.join(sys.argv[3:])+"nZ"+str(nZ)+"nN"+str(nN))
