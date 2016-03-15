import numpy as np
import matplotlib.pyplot as plt
import sys

def drawPsi(psiActual, xGrid, zGrid, name, graphtype, clearAfter = True):
    # graphtype = 'streamplot' gives a streamplot of the fluid velocity derived from psi
    #           = 'contour' gives a contour map of psi
    
    validGraphTypes = ['streamplot', 'contour']

    if graphtype == 'contour':
        try:
            plt.contour(xGrid, zGrid, psiActual, 6, colors='k')
        except ValueError:
            pass
    if clearAfter:
        plt.cla()

def drawTemp(tempActual, xGrid, zGrid, name, graphType,              xLabel = '', zLabel = '', clearAfter = True):
    # Possible graph types:
    # graphType == 'heatmap' draws a heatmap
    # == 'contour' draws a contour map
    
    validGraphTypes = ['heatmap', 'contour']
    
    if xLabel:
        plt.xlabel(xLabel)
    if zLabel:
        plt.ylabel(zLabel)
#    plt.axes().set_aspect('equal')

    if graphType == 'heatmap':
        plt.pcolormesh(xGrid, zGrid, tempActual, cmap='coolwarm')
    elif graphType == 'contour':
        plt.contour(xGrid, zGrid, tempActual)
    
    if graphType in validGraphTypes:
        pass
    else:
        print "Invalid graph type: " + graphType
        print "Should be one of: " + ", ".join(validGraphTypes)
    if clearAfter:
        plt.cla()


def main():
    a = 3
    Nn = 101
    Nz = 101
    Nx = int(a*Nz)
    xAxis = np.linspace(0, a, num=Nx)
    zAxis = np.linspace(0, 1, num=Nz)
    xGrid, zGrid = np.meshgrid(xAxis,zAxis)
    data = np.fromfile(sys.argv[1], dtype=np.dtype(np.double))
    temp = np.transpose(data[:Nn*Nz].reshape(Nn, Nz))
    psi = np.transpose(data[2*Nn*Nz:3*Nn*Nz].reshape(Nn, Nz))
    #plt.gca().set_aspect(1.0)
    #plt.axis('off')
    #plt.gca().get_xaxis().set_visible(False)
    #plt.gca().get_yaxis().set_visible(False)
    #plt.plot(temp[:, 1])

    tempActual = np.zeros((Nz, Nx))
    psiActual = np.zeros((Nz, Nx))
    #for n in xrange(Nn):
    #    for k in xrange(Nz):
    #        tempActual[k, :] += temp[k, n]*np.cos(n*np.pi*xAxis)
    preCos = np.cos(np.pi/a*np.outer(np.arange(Nn), xAxis))
    preSin = np.sin(np.pi/a*np.outer(np.arange(Nn), xAxis))
    np.dot(psi, preSin, out=psiActual)
    np.dot(temp, preCos, out=tempActual)
    
    plt.subplot(1, 2, 1)
    plt.gca().set_aspect(1.0)
    plt.yticks([0, 1])
    plt.xticks([0, 1, 2, 3])
    plt.tick_params(axis='x', bottom='off', top='off')
    drawTemp(tempActual, xGrid, zGrid, sys.argv[1], 'heatmap', clearAfter=False)
    drawTemp(tempActual, xGrid, zGrid, sys.argv[1], 'contour', clearAfter=False)
    #plt.savefig(sys.argv[1] +'tmp.png', bbox_inches='tight', pad_inches=0)
    plt.subplot(1, 2, 2)
    plt.gca().set_aspect(1.0)
    plt.gca().get_yaxis().set_visible(False)
    plt.xticks([0, 1, 2, 3])
    plt.tick_params(axis='x', bottom='off', top='off')
    plt.gcf().tight_layout()
    drawPsi(psiActual, xGrid, zGrid, sys.argv[1], 'contour', clearAfter=False)
    plt.savefig(sys.argv[1] +'.png', bbox_inches='tight', pad_inches=0, dpi=200)
    #plt.show()

main()
