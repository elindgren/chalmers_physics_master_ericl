# Internal imports
import os

# External imports
import numpy as np
import matplotlib.pyplot as plt


# Set plot params
plt.rc('font', size=14)          # controls default text sizes
plt.rc('axes', titlesize=14)     # fontsize of the axes title
plt.rc('axes', labelsize=14)     # fontsize of the x and y labels
plt.rc('xtick', labelsize=14)    # fontsize of the tick labels
plt.rc('ytick', labelsize=14)    # fontsize of the tick labels
plt.rc('legend', fontsize=14)    # legend fontsize


'''
    This script calculates the linear expansion coefficient in Fe
    from MD simulations at various temperatures.

    Devloped by Eric Lindgren, SNUID: 2020-81634
    June 2020
'''
l0 = 2.8553 # Å
lx = [l0]
T = [0]

# Extract final lattice constant for each temperature
for fl in os.listdir('.'):
    if 'lx' in fl:
        ts = fl.split('_')[1]
        with open(fl, 'r') as f:
            r = f.readlines()[-1].rstrip()
            l = float(r.split(' ')[1])
            lx.append(l)
        # Get corresponding temperature
        for ft in os.listdir('.'):
            if f'T_{ts}' in ft:
                with open(ft, 'r') as d:
                    r = d.readlines()[-1].rstrip()
                    t = float(r.split(' ')[1])
                    T.append(t)
lx = np.array( lx )
T = np.array( T )
sort = np.argsort(lx)
lx = lx[sort]
T = T[sort]

# Perform linefit
norm_lx = lx - l0 # normalize
norm_lx /= l0

fit = np.polyfit( x=T, y=norm_lx, deg=2 )
lexpc = 2*fit[0] # 1/K^2 - multiply by two since we want the derivative dl/dT
# Calculate elastic constant
print(f'Linear expansion coefficient: {lexpc:e} 1/K^2')

# Plot stress-strain
fig, ax = plt.subplots(figsize=(8,6))
ax.plot( T, norm_lx, linewidth=2, label='Data' )
Tp = np.linspace(T[0], T[-1], 100)
ax.plot(Tp, fit[0]*Tp**2 + fit[1]*Tp + fit[2], label='Fit')
ax.set_xlabel(r'$T$, K')
ax.set_ylabel(r'$l_x/l_0$, (Å)')
ax.grid()
ax.legend(loc='best')
plt.tight_layout()
plt.savefig('lin_expansion.png')
plt.show()
