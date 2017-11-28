import numpy as np
import sys

nZ = int(sys.argv[1])
nN = int(sys.argv[2])
n = int(sys.argv[3])
amp = float(sys.argv[4])

inFile = sys.argv[5]
outFile = sys.argv[6]

print "Adding " + str(amp) + " perturbation to mode " + str(n) + "; " + inFile + " -> " + outFile

data = np.fromfile(inFile, dtype=np.dtype(np.double))

data[n*nZ:(n+1)*nZ] += amp*np.sin(np.pi*np.linspace(0, 1, nZ))

data.tofile(outFile)

