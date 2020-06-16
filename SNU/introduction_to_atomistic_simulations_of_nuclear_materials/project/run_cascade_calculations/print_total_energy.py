import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from cycler import cycler

# Set plot params
plt.rc('font', size=16)          # controls default text sizes
plt.rc('axes', titlesize=16)     # fontsize of the axes title
plt.rc('axes', labelsize=16)     # fontsize of the x and y labels
plt.rc('xtick', labelsize=16)    # fontsize of the tick labels
plt.rc('ytick', labelsize=16)    # fontsize of the tick labels
plt.rc('legend', fontsize=14)    # legend fontsize
line_cycler = cycler('linestyle',['-','--',':','-.', '-', '--']) + cycler('color', ['b', 'r', 'g', 'k', 'c', 'y'])
plt.rc('axes', prop_cycle=line_cycler)


df = pd.read_csv('runstats.out', skiprows=1, names=['time','etotal','press','temp','c_kemax','c_emax'])
etot = df['etotal']

fig, ax = plt.subplots(figsize=(8,6))
ax.plot(df['time'], etot-etot[0], linewidth=2)
ax.set_xlabel('Time (ps)')
ax.set_ylabel('Change in total energy (eV)')
ax.grid()
plt.savefig('total_energy.png')