# Internal imports
import os 
from cycler import cycler

# External imports
import numpy as np
import pandas as pd
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

E = []
bulk = {
    'mu' : [],
    'sigma' : [] 
}
vol = {
    'mu' : [],
    'sigma' : [] 
}

for obj in os.listdir('.'):
    if 'eV' in obj:
        E.append(float(obj.split('_')[0]))
        df = pd.read_csv(obj, sep=',', names=['direction','bulk','vol'])
        # Bulk modulus
        bulk['mu'].append(np.mean(df['bulk']))
        bulk['sigma'].append(np.std(df['bulk']))
        # Volume
        vol['mu'].append(np.mean(df['vol']))
        vol['sigma'].append(np.std(df['vol']))

# Scatter results
print(bulk)
fig, ax= plt.subplots(figsize=(8,6))
ax.errorbar(x=E, y=bulk['mu'], yerr=bulk['sigma'], fmt='o', ms=10, capsize=10, capthick=2, label='Bulk modulus')
ax.set_xlabel('PKA energy (eV)')
ax.set_ylabel('Bulk modulus (GPa)')
ax.grid()
ax.legend(loc='best')
plt.tight_layout()
plt.show()

