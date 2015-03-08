import numpy as np
import matplotlib.pyplot as plt

with open('/home/manyu/BTP/cuda_code/file_xelec.dat') as file:
    z = [[float(digit) for digit in line.split()] for line in file]

fig, axs = plt.subplots()
z1 = z[50]
x = np.linspace(1, 100, 100)
X, Y = np.meshgrid(x, x)

maxz = max(max(z))
minz = min(min(z))

levels = np.linspace(minz-1, maxz+1, 200)

#cs = axs[0].contourf(X, Y, z, levels=levels)
#fig.colorbar(cs)

axs.plot(x,z1)

plt.title('1D Electric Field Plot')
plt.xlabel('x grid dimension')
plt.ylabel('Electric Field (V/m)')


plt.show()