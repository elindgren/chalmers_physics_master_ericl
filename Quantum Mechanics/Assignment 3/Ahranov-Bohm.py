import numpy as np
import matplotlib.pyplot as plt

def energy_eig(l, j):
    # Define j = q*B0*rho^2/2hbar
    return (l-j)**2

# Constant in front of E: c = hbar^2/(2*m*rho^2)

ls = np.arange(-5,5)
j = np.linspace(-10,10,100)

fig, ax = plt.subplots()
plasma = plt.cm.plasma(np.linspace(0, 1, 10))

for l in ls:
    ax.plot(j, energy_eig(l,j), label=f'l={l}', color=plasma[l])

ax.set_xlabel("q*B0*rho^2/2hbar")
ax.set_ylabel("proportional to E_l")
plt.title("Plot of energy eigenvalues E_l as a function of enclosed flux")
plt.grid()
plt.legend()
plt.savefig("ahranov_bohm.png")
plt.show()