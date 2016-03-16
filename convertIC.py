import numpy as np
import sys
nZOld = sys.argv[1]
nNOld = sys.argv[2]
nZNew = sys.argv[3]
nNNew = sys.argv[4]
infile = sys.argv[5]
outfile = sys.argv[6]

zRatio = float(nZNew - 1)/(nZOld - 1)

dataOld = np.fromfile(infile, dtype=np.dtype(np.double)
dataNew = np.zeros(nZNew*nNNew)

for n in xrange(nNOld):
	# interpolate end points
	dataNew[n*nZNew+0] = dataOld[n*nZOld+0] 
	dataNew[n*nZNew+nZNew] = dataOld[n*nZOld+nZOld] 

dataNew.tofile(outfile)

