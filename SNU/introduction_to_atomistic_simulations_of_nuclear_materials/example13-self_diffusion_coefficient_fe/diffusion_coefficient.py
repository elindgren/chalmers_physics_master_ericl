# Internal imports
from cycler import cycler

# External imports
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# Set plot params
plt.rc('font', size=16)          # controls default text sizes
plt.rc('axes', titlesize=16)     # fontsize of the axes title
plt.rc('axes', labelsize=16)     # fontsize of the x and y labels
plt.rc('xtick', labelsize=16)    # fontsize of the tick labels
plt.rc('ytick', labelsize=16)    # fontsize of the tick labels
plt.rc('legend', fontsize=14)    # legend fontsize
line_cycler = cycler('linestyle',['-','--',':','-.', '-', '--']) + cycler('color', ['b', 'r', 'g', 'k', 'c', 'y'])
plt.rc('axes', prop_cycle=line_cycler)


cases = ['perfect', 'SIA'] # Define cases
d = 3 # Dimensionality

fig, ax = plt.subplots(figsize=(8,6))
for case in cases:
    df = pd.read_csv(f'msd_sc4-{case}.out', skiprows=1, names=['t', 'msdx', 'msdy', 'msdz', 'msd'])
    msd = np.array(df['msd'])
    t = np.array(df['t'])

    # Calculate self diffusion coefficient
    fit = np.polyfit(t, msd, deg=1)
    D = fit[0]/(2*d) # Å^2/ps
    print(f'Case={case}: D = {D:e} Å^2/ps')

    ax.plot(t, msd, linewidth=2, label=f'Case: {case}')
ax.set_xlabel('Time, (ps)')
ax.set_ylabel(r'$\bar{r^2}$, ($\rm Å^2$)')
ax.grid()
ax.legend(loc='best')
plt.show()
