import numpy as np
import matplotlib.pyplot as plt

# Set plot params
plt.rc('font', size=14)          # controls default text sizes
plt.rc('axes', titlesize=14)     # fontsize of the axes title
plt.rc('axes', labelsize=14)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=14)    # fontsize of the tick labels
plt.rc('ytick', labelsize=14)    # fontsize of the tick labels
plt.rc('legend', fontsize=14)    # legend fontsize

def free_energy(m, T_tilde, Tc, kB, J, z):
    T = Tc*(T_tilde+1)
    return -0.5*J*z*m**2 + 0.5*kB*T*( (1+m)*np.log(1+m) + (1-m)*np.log(1-m) - np.log(2))  # It was kinda fine until I took -np.log(2)


# Constants
kB = 1
J = 1
z = 4  # 2D lattice
Tc = z*J/kB

T_tilde = [-1, -0.5, 0, 0.5, 1]
style = ['-', '--', '-.', ':', '-']
m = np.linspace(-5,5,99)

fig, ax = plt.subplots(figsize=(10,6))

for i, Tt in enumerate(T_tilde):
    F = free_energy(m, Tt, Tc, kB, J, z)
    s = style[i]
    ax.plot(m, F, linestyle=s, c='C'+str(2*i), linewidth=3, label=r'$\hat{T}$ = ' + f'{Tt}')

ax.set_xlabel(r'Magnetization $m$ (T)')
ax.set_ylabel(r'Free energy, $F$ (J)')
ax.grid()

ax.legend(loc='best')
plt.show()
