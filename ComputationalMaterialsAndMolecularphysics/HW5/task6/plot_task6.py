# Internal imports
import pickle

# External imports
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

# ASE
from ase import Atoms
from ase.db import connect


# Set plot params
plt.rc('font', size=18)          # controls default text sizes
plt.rc('axes', titlesize=18)     # fontsize of the axes title
plt.rc('axes', labelsize=18)     # fontsize of the x and y labels
plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
plt.rc('ytick', labelsize=18)    # fontsize of the tick labels
plt.rc('legend', fontsize=18)    # legend fontsize

# Connect to DB
bulkDB = connect('./bulk.db')

# Extract and plot convergence data
fig, ax = plt.subplots(figsize=(8,6))
dbList = list(bulkDB.select())

ks = dbList[0].data['ks']
E = dbList[0].data['energies']
ax.plot(ks, E, linewidth=2)

ax.set_xlabel(r'Number of $k$-points')
ax.set_ylabel('Energy (eV)')
ax.grid()
plt.tight_layout()
plt.savefig('kConvergence.png')

# Extract and plot band electronic band structure and DOS
bs = pickle.load(open( "Ebs.p", "rb" ))
d = pickle.load(open( "Edos.p", "rb" ))

fig, ax = plt.subplots(1,2, figsize=(12,6))

# BS
# bs.energies = bs.energies - d['fermi']
emax = 16
emin = -13
bs.plot(filename='', ax=ax[0], show=False, emax=emax, emin=emin)
ax[0].set_ylabel(r'Energy (eV)')

# DOS
# e = d['e']-d['fermi']
e = d['e']
# ax[1].fill_between( d['dos'], e, y2=0, color='grey',
#                    edgecolor='k', lw=1, alpha=0.6, label=r'$g(E)$')
ax[1].plot(d['dos'], e, label=r'$g(E)$' )

ePos = np.array([ ei for ei in e if ei>=0 ])
# Calculate free electron density
# Na = 6.02214076e23  # Avogadro's constant
# Z = 3 # Nbr of valence electrons of Al
# rho = 2720 # kg/m3 
# ma = 26.98 * 1.66e-27
# n = Na*Z*rho/ma
# freeE = 1.5 * n/d['fermi'] * np.sqrt(ePos/d['fermi'])   
freeE = 0.061*np.sqrt(ePos)  # TODO set proper scale
ax[1].plot(freeE, ePos-11.2, color='k', label=r'Free $\rm e^{-}$, $g(E) \propto \sqrt{E}$')

ax[1].legend(loc='best')
ax[1].set_xticks([])
# ax[1].set_xlabel(r'DOS ($\rm eV^{-1}$)') # TODO set proper units
ax[1].set_xlabel(r'DOS')
ax[1].set_ylabel(r'Energy (eV)')
# ax[1].set_ylim((emin, emax))
plt.tight_layout()
plt.savefig('electronicAl.png')
