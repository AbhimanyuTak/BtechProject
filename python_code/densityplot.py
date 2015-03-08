import numpy as np
import matplotlib.pyplot as plt

with open('/home/manyu/BTP/cuda_code/file_a.dat') as file:
    z = [[float(digit) for digit in line.split()] for line in file]
    
maxz = max(max(z))
minz = min(min(z))

fig, axs = plt.subplots()

x = np.linspace(1, 100, 100)
X, Y = np.meshgrid(x, x)

levels = np.linspace(minz, maxz, 5)

cs = axs.contourf(X, Y, z, levels=levels)
cbar = fig.colorbar(cs)
#cbar.set_label('Density')
plt.title('Density Plot')
plt.xlabel('x grid dimension')
plt.ylabel('y grid dimension')

plt.show()