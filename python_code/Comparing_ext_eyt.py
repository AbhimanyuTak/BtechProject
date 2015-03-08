import matplotlib.pyplot as plt
import numpy as np

x = []
y = []
z = []
c = 0

# Simple data to display in various forms
#x = np.linspace(0, 2 * np.pi, 400)
#y = np.sin(x ** 2)

readFile = open('/home/manyu/BTP/cuda_code/eyt.dat','r')
sepFile = readFile.read().split('\n')
readFile.close()

for plotPair in sepFile:
    xAndY = plotPair.split(' ')
    c=c+1
    try:
    	x.append(float(xAndY[0]))
        y.append(float(xAndY[1]))
    	z.append(float(xAndY[2]))
    except ValueError,e:
        print "error",e,"on line",c

plt.close('all')

n = 1000
# Two subplots, the axes array is 1-d
f, axarr = plt.subplots(2, sharex=True)
axarr[0].plot(x[1:n], y[1:n])
axarr[0].set_title('Electric Field in X Direction after Scattering')
axarr[1].plot(x[1:n], z[1:n])
axarr[1].set_title('Electric Field in Y Direction')

plt.xlabel('Time Steps')
plt.ylabel('Electric Field (V/m)')

plt.show()