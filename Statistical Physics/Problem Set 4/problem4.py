import numpy as np
import matplotlib.pyplot as plt
import scipy.special as spec
plt.rc('font', size=18) # controls default text sizes
plt.rc('axes', titlesize=18) # fontsize of the axes title
plt.rc('axes', labelsize=18) # fontsize of the x and y labels
plt.rc('xtick', labelsize=18) # fontsize of the tick labels
plt.rc('ytick', labelsize=18) # fontsize of the tick labels
plt.rc('legend', fontsize=18) # legend fontsize


def C(x,t, b, c, eta):
    return ( 2*np.cosh(np.sqrt(b/c)*x) - np.exp(np.sqrt(b/c)*x)*spec.erf(np.sqrt(b*eta*t) + x/np.sqrt(4*eta*c*t)) - np.exp(-np.sqrt(b/c)*x)*spec.erf(np.sqrt(b*eta*t) - x/np.sqrt(4*eta*c*t)) )

# Constants
b = 1
c = 1
eta = 1
ts = [1e-5, 1e-1]
ls = ['--', ':']
L = 1
x = np.linspace(-L,L,100)

fig, ax = plt.subplots(figsize=(8,6))

for i, t in enumerate(ts): 
    func = C(x, t, b, c, eta)
    ax.plot(x, func/np.max(func), linewidth=3, linestyle=ls[i], label=rf'$t$={t:.2E}')
ax.grid()
ax.legend(loc='best')
ax.set_xlabel(r'$x$, arb. units')
ax.set_ylabel(r'$C(x,t)$, arb. units')
plt.savefig("problem4.png")
plt.show()
