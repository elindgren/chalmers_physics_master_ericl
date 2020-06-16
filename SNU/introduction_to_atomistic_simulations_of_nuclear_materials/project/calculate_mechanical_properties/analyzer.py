# Internal imports
import os 
from cycler import cycler

# External imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# Set plot params
plt.rc('font', size=16)          # controls default text sizes
plt.rc('axes', titlesize=20)     # fontsize of the axes title
plt.rc('axes', labelsize=24)     # fontsize of the x and y labels
plt.rc('xtick', labelsize=20)    # fontsize of the tick labels
plt.rc('ytick', labelsize=20)    # fontsize of the tick labels
plt.rc('legend', fontsize=18)    # legend fontsize
line_cycler = cycler('linestyle',['-','--',':','-.', '-', '--']) + cycler('color', ['b', 'r', 'g', 'k', 'c', 'y'])
plt.rc('axes', prop_cycle=line_cycler)

E = []
mech = {
    'bulk' : {
        'perfect' : 310.02871485826,
        'label' : 'Bulk modulus (GPa)',
        'mean' : [],
        'std' : []
    },
    'shear2' : {
        'perfect' : 159.050933744743,
        'label' : r'$\rm (C_{11}+C_{12})/2$ (GPa)',
        'mean' : [],
        'std' : []
    },
    'shear1' : {
        'perfect' : 161.063454093376,
        'label' : r'$\rmC_{44}$ (GPa)',
        'mean' : [],
        'std' : []
    },
    'poisson' : {
        'perfect' : 0.2809491519789,
        'label' : 'Poisson ratio',
        'mean' : [],
        'std' : []
    },
}

for obj in os.listdir('.'):
    if 'eV' in obj:
        E.append(float(obj.split('_')[0]))
        df = pd.read_csv(obj, sep=',', names=['direction','bulk','shear1', 'shear2', 'poisson'])
        for col, data in df.iteritems():
            if col != 'direction':
                mech[col]['mean'].append(np.mean(data))
                mech[col]['std'].append(np.std(data))

# Scatter results
for prop, data in mech.items():
    fig, ax= plt.subplots(figsize=(12,9))
    ax.errorbar(x=E, y=data['mean'], yerr=data['std'], fmt='o', ms=10, capsize=10, capthick=2, label='Defect system')
    ax.axhline(data['perfect'], c='k', linestyle='--', linewidth=2, label='Perfect system')
    ax.set_xlabel('PKA energy (eV)')
    ax.set_ylabel(data['label'])
    ax.grid()
    ax.legend(loc='best')
    plt.tight_layout()
    plt.savefig(f'figures/{prop}.png')

