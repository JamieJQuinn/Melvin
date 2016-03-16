import numpy as np
import sys
nZOld = int(sys.argv[1])
nNOld = int(sys.argv[2])
nZNew = int(sys.argv[3])
nNNew = int(sys.argv[4])
infile = sys.argv[5]
outfile = sys.argv[6]

if (nZOld > nZNew) or (nNOld > nNNew):
    print "Converting to a lower resolution. Not allowed."
    exit(0)

print "Converting " + infile + " (" + str(nZOld) + ", " + str(nNOld) + ") to " + outfile + " (" + str(nZNew) + ", " + str(nNNew) + ")"

zRatio = float(nZNew - 1)/(nZOld - 1)

dataOld = np.fromfile(infile, dtype=np.dtype(np.double))
dataNew = np.zeros(3*nZNew*nNNew)

for i in xrange(3):
    for n in xrange(nNOld):
            # interpolate 
            dataNew[(i*nNNew+n)*nZNew:(i*nNNew+n+1)*nZNew] = np.interp(np.linspace(0, 1, num=nZNew), np.linspace(0, 1, num=nZOld), dataOld[(i*nNOld+n)*nZOld:(i*nNOld+n+1)*nZOld])


dataNew.tofile(outfile)

