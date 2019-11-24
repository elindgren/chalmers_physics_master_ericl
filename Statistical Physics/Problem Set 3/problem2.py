import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns



# Set plot params
plt.rc('font', size=14)          # controls default text sizes
plt.rc('axes', titlesize=14)     # fontsize of the axes title
plt.rc('axes', labelsize=14)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=14)    # fontsize of the tick labels
plt.rc('ytick', labelsize=14)    # fontsize of the tick labels
plt.rc('legend', fontsize=14)    # legend fontsize


def free_energy(nA, b):
    ''' Calculate the value of the free energy for a given nA and b(eta).'''
    return -0.5*( nA**2 + 0.5*(1-nA)**2 ) + b*( nA*np.log(nA) + (1-nA)*np.log(0.5 * (1-nA)) )


# Define plot range
nA = np.linspace(0.00001,0.99999,100)
b = np.linspace(0,100,100)

N, B = np.meshgrid(nA, b)
F = free_energy(N, B)

# Define plots
fig = plt.figure(figsize=(10,6))
ax = fig.add_subplot(111, projection='3d')
ax.contour3D(N, B, F, 50, cmap='binary')
# ax.plot3D(nA, b, free_energy(nA, b))

plt.show()


