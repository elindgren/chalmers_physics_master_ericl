import numpy as np
import matplotlib.pyplot as plt

def energy_eig(l, j):
    # Define j = q*B0*rho^2/2hbar
    return (l-j)**2

# Constant in front of E: c = hbar^2/(2*m*rho^2)

ls = np.arange(0,10)
j = np.linspace(0,10,100)

fig, ax = plt.subplots()

for l in ls:
    ax.plot(j, energy_eig(l,j), label=f'l={l}')

plt.legend()
plt.show()