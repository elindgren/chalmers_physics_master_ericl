# Internal imports
from cycler import cycler

# External imports
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


# Set plot params
plt.rc('font', size=16)          # controls default text sizes
plt.rc('axes', titlesize=16)     # fontsize of the axes title
plt.rc('axes', labelsize=16)     # fontsize of the x and y labels
plt.rc('xtick', labelsize=16)    # fontsize of the tick labels
plt.rc('ytick', labelsize=16)    # fontsize of the tick labels
plt.rc('legend', fontsize=14)    # legend fontsize
line_cycler = cycler('linestyle',['-','--',':','-.', '-', '--']) + cycler('color', ['r', 'g', 'b', 'k', 'c', 'y'])
plt.rc('axes', prop_cycle=line_cycler)


'''
    This script integrates the VACF to calcualte MSD according to the Green-Kubo relations.
'''

# Load data
vacdf = pd.read_csv('vacf.out', skiprows=1, names=['t', 'x', 'y', 'z', 'tot','msd_3d'])

# Calculate MSD - just trapz
t = vacdf['t']
msdx = np.zeros(len(t))
msdy = np.zeros(len(t))
msdz = np.zeros(len(t))
msdtot = np.zeros(len(t))
for i in range(len(t)):
    msdx[i] = np.trapz(y=vacdf['x'][:i], x=t[:i])
    msdy[i] =  np.trapz(y=vacdf['y'][:i], x=t[:i])
    msdz[i] =  np.trapz(y=vacdf['z'][:i], x=t[:i])
    msdtot[i] =  np.trapz(y=vacdf['tot'][:i], x=t[:i]) # Å^2/ps

fit = np.polyfit(t, vacdf['msd_3d'], deg=1)
print(f'Reference diffusion coefficient: {fit[0]/6:e} Å^2/ps')
idx = np.argwhere(t==2.0)[0]
print(f'Diffusion coefficient at 2 ps: {msdtot[idx[0]]:e} Å^2/ps')

# Plot
fig, ax = plt.subplots(figsize=(8,6))
ax.plot(t, msdtot, linewidth=2, label=r'$\rm D{tot}$')
ax.plot(t, msdx, linewidth=2, label=r'$\rm D_{x}$')
ax.plot(t, msdy, linewidth=2, label=r'$\rm D_{y}$')
ax.plot(t, msdz, linewidth=2, label=r'$\rm D_{z}$')
ax.set_xlabel('Time, (ps)')
ax.set_ylabel('D, (Å^2/ps)')
ax.legend(loc='best')
ax.set_xlim(0,2) # Only the first 2 ps are good before errors starts to accumulate too much - pretty converged at 1 ps
ax.grid()
plt.tight_layout()
plt.show()


