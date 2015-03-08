import numpy as np
import matplotlib.pyplot as plt

with open('/home/manyu/BTP/cuda_code/file_xelec.dat') as file:
    z = [[float(digit) for digit in line.split()] for line in file]

fig, axs = plt.subplots()

x = np.linspace(1, 100, 100)
X, Y = np.meshgrid(x, x)

#maxz = max(max(z))
#minz = min(min(z))
maxz = np.amax(z)
minz = np.amin(z)

levels = np.linspace(minz-1, maxz, 150)

cs = axs.contourf(X, Y, z, levels=levels)
cbar = fig.colorbar(cs)
cbar.set_label('V/m')
plt.title('Electric Field')
plt.xlabel('x grid dimension')
plt.ylabel('y grid dimension')

plt.show()