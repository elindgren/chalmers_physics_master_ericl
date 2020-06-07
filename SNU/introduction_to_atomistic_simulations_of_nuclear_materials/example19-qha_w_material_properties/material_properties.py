# Internal imports
from cycler import cycler

# External imports
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


df = pd.read_table('w.output', sep='    ', skiprows=2, header=None, names=['T', 'E_0', 'B_0', "B'_0", 'V0'])


# Plot bulk modulus as a function of temperature
fig, ax = plt.subplots(figsize=(8,6))
ax.plot(df['T'], df['B_0'], linewidth=2)
ax.grid()
ax.set_xlabel('Temperature, (K)')
ax.set_ylabel('Bulk modulus, (GPa)')
plt.show()