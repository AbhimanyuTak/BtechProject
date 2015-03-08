import numpy as np
from numpy import arange, sqrt, sin
from tvtk.api import tvtk
from mayavi.scripts import mayavi2
import matplotlib.pyplot as plt
from mayavi.mlab import *
from mayavi import mlab

readFile = open('eyt.dat','r')
sepFile = readFile.read().split('\n')
readFile.close()

#readFile = open('eyi1.dat','r')
#sepFile1 = readFile.read().split('\n')
#readFile.close()

#x, y = np.loadtxt(sepFile, delimiter=' ', usecols=(0, 1), unpack=True)
#z = x/x*0.00000001

x, y, z= np.loadtxt(sepFile, delimiter=' ', usecols=(0, 1, 2), unpack=True)
#z = x/x*0.01
#z = np.sqrt(x**2 + y**2)

n = 10000
s = plot3d(x[1:n], y[1:n], color=(1,1,1), colormap='Spectral', tube_radius=0.5)
#r = plot3d(a[1:n], b[1:n], c[1:n], color=(1,0,1), colormap='Spectral', tube_radius=1.0)

s.show()
#r.show()

