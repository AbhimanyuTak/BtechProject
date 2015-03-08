from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np

u = []
v = []
w = []
x = []
y = []
z = []
c = 0

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

#for i in u:
#    u[i] = np.sqrt(y[i]**2 + z[i]**2)
a = 1
b = 2000
plt.plot(x[a:b],z[a:b])
plt.title('Electric Field in Y direction')
plt.xlabel('x axis label')
plt.ylabel('y axis label')
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#X, Y, Z = axes3d.get_test_data(0.05)
#ax.plot_wireframe(y[a:b], z[a:b], u[a:b], rstride=10, cstride=10)

plt.show()