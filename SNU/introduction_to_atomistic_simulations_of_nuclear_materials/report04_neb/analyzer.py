# Internal imports
from cycler import cycler

# External imports
import numpy as np
import matplotlib.pyplot as plt

# Set plot params
plt.rc('font', size=16)          # controls default text sizes
plt.rc('axes', titlesize=16)     # fontsize of the axes title
plt.rc('axes', labelsize=16)     # fontsize of the x and y labels
plt.rc('xtick', labelsize=16)    # fontsize of the tick labels
plt.rc('ytick', labelsize=16)    # fontsize of the tick labels
plt.rc('legend', fontsize=14)    # legend fontsize
linestyle_cycler = cycler('linestyle',['-','--',':','-.'])
plt.rc('axes', prop_cycle=linestyle_cycler)

'''
    This script takes the output from each case in runner.py and 
    plots the transition energy as a function of reaction coordinate
    and calculates the transition barrier. 

    Developed by Eric Lindgren, SNU ID ericlin
    May 2020
'''

trans = {}

# Extract transition details for each case
with open('transitions.txt', 'r') as f:
    for i, line in enumerate(f):
        if 'case' in line and i == 0:
            c = line.rstrip()
            t = []
        elif 'case' in line and i != 0:
            # Save old data structures
            t = np.array( t )
            trans[c] = t
            # Create new ones
            c = line
            t = []
        else:
            sl = line.split(' ')
            r = float(sl[0])
            e = float(sl[1])
            t.append([r, e])
    t = np.array( t )
    trans[c] = t

# Plot transition
fig, ax = plt.subplots(figsize=(8,6))
for key, val in trans.items():
    r = val[:,0]
    n_e = val[:,1] - val[0,1]
    assert val[0,1] == val[-1,1] # These should match for a good NEB calculation
    lab = f'Case ({key[-1]})'
    ax.plot(r, n_e, markersize=2, linewidth=2, label=lab)
    # Find energy barrier
    print(f'{lab}: Em={n_e.max():.3f} eV')
ax.legend(loc='best')
ax.set_ylabel(r'Transition energy $E_m$ (eV)')
ax.set_xlabel(r'Reaction coordinate $R$')
ax.grid()
plt.tight_layout()
plt.savefig('SIA_transitions.png')