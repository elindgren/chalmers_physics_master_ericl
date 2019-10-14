import numpy as np
import matplotlib.pyplot as plt

def energy_eig(l, js):
    # Define j = q*B0*rho^2/2hbar
    E_vec = np.zeros(len(js))
    for idx, j in enumerate(js):
        #if not j < l:
        E_vec[idx] = (l-j)**2

    return E_vec

# Constant in front of E: c = hbar^2/(2*m*rho^2)

ls = np.arange(-4,4)
j = np.linspace(-10,10,100)

fig, ax = plt.subplots()

for l in ls:
    ax.plot(j, energy_eig(l,j), label=f'l={l}')

plt.grid()
plt.legend()
plt.show()