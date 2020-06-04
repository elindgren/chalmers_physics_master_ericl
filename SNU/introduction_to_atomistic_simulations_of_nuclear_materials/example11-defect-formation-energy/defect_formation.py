# Internal imports
import pickle
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
line_cycler = cycler('linestyle',['-','--',':','-.', '-', '--']) + cycler('color', ['b', 'r', 'g', 'k', 'c', 'y'])
plt.rc('axes', prop_cycle=line_cycler)

# Load Ed data
with open('Ed.pckl', 'rb') as f:
    Ed = pickle.load(f)


fig, ax = plt.subplots(3,1, figsize=(10,10))
# Calculate formation energies:
for s, data in Ed.items():
    e_perf = data['perfect']['E']
    n_perf = data['perfect']['N']

    e_v1 = data['V1']['E']
    n_v1 = data['V1']['N']

    e_v2 = data['V2']['E']
    n_v2 = data['V2']['N']
    
    e_sia = data['SIA']['E']
    n_sia = data['SIA']['N']

    E_v1 = n_v1 * ( e_v1/n_v1 - e_perf/n_perf)
    E_v2 = n_v2 * ( e_v2/n_v2 - e_perf/n_perf)
    E_sia = n_sia * ( e_sia/n_sia - e_perf/n_perf)

    ax[0].scatter(s, E_v1)
    ax[1].scatter(s, E_v2)
    ax[2].scatter(s, E_sia)    

ax[0].set_xlabel('Supercell size')
ax[0].set_ylabel('V1 Defect formation energy, (eV)')
ax[1].set_xlabel('Supercell size')
ax[1].set_ylabel('V2 Defect formation energy, (eV)')
ax[2].set_xlabel('Supercell size')
ax[2].set_ylabel('SIA Defect formation energy, (eV)')
plt.tight_layout()
plt.show()
