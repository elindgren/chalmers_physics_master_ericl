import numpy as np
import matplotlib.pyplot as plt
plt.rc('font', size=18)          # controls default text sizes
plt.rc('axes', titlesize=18)     # fontsize of the axes title
plt.rc('axes', labelsize=18)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
plt.rc('ytick', labelsize=18)    # fontsize of the tick labels
plt.rc('legend', fontsize=18)    # legend fontsize

# Define function

def C(l, t, A, L, R, c):
    k = 2*np.pi*l/L
    return l*(l+1)*A/k**2 * np.abs( np.cos( k * c * np.sqrt(1/(3*(1+R))) * t ) )**2 

# Constants
c = 1  # ly/y
t = 380000  # 380000 years
A = 1
# L = 9.7  # ly
R = 0.7
L = 440 * c * t / np.sqrt((3*(1+R)))
print(f'L={L:.4f} ly')
ls = np.linspace(40, 1500, 1000)

# plot
fig, ax = plt.subplots(figsize=(8,6))
ax.plot(ls, C(ls, t, A, L, R, c))
ax.set_xlabel(r'$\ell$')
ax.set_ylabel(r'$\hat{C}(x,y)$, (ly$^2$)')
plt.savefig("problem1.png")
plt.show()