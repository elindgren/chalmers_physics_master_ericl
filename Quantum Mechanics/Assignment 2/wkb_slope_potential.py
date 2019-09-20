import numpy as np
import matplotlib.pyplot as plt
from decimal import Decimal


def energy_equation(energy, n):
    tmp = np.sqrt(2 * m) / hbar * (-L / V0) * 2 / 3 * (energy ** (2 / 3) - (energy + V0) ** (2 / 3)) - (np.pi/4+np.pi*n)
    #print(tmp)
    return tmp

# Define constants
hbar = 1
m = 1
# Parameters
V0 = 1  # 100 eV
L = 1

# Define energy space
E = np.linspace(0, 100, 100)

# Target points: pi/4 + n*pi
n = np.array(range(20))  # 20 target points
target = np.pi/4*np.ones((len(n),)) + np.pi*n

# Define function vector
fun = np.zeros((len(target), 1))

# Define plot
fig, ax = plt.subplots()

# Functions which are to be zerod
for i in range(len(n)):
    equation_n = np.zeros((len(E), 1))
    for j in range(len(E)):
        equation_n[j] = energy_equation(E[j], n[i])
    # For each n, calculate fun - n*pi/2+pi/4, and find it's zero
    # Solve equation
    # TODO
    # Plot function
    ax.plot(E, equation_n, 'b')
    # Plot solution to equation

plt.xlabel("Energy, E")
plt.ylabel("fun-n*pi/2+pi/4")
plt.show()

# TODO reduce to dimensionless constants
