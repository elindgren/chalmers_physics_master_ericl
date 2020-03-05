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
ax.plot(ks, E)

ax.set_xlabel(r'Number of $k$-points')
ax.set_ylabel('Energy (eV)')
ax.grid()
plt.tight_layout()
plt.savefig('convergenceSi')

# Extract and plot band electronic band structure and DOS
bs = pickle.load(open( "Ebs.p", "rb" ))
d = pickle.load(open( "Edos.p", "rb" ))

fig, ax = plt.subplots(1,2, figsize=(12,6))
# BS
bs.plot(filename='Al_Electronic_bs.png', ax=ax[0], show=False, emax=10.0)
# DOS
ax[1].plot(d['e']-d['fermi'], d['dos'])
ax[1].set_xlabel(r'Energy relative to $\epsilon_F$ (eV)')
ax[1].set_ylabel('DOS') # TODO set proper units
plt.savefig('electronicSi')

# Extract and plot phonon band structure and DOS
bs = pickle.load(open( "Pbs.p", "rb" ))
dos = pickle.load(open( "Pdos.p", "rb" ))

fig, ax = plt.subplots(1,2, figsize=(12,6))
emax = 0.05
bs.plot(ax=ax[0], emin=0.0, emax=emax)

ax[1].fill_between(dos.weights[0], dos.energy, y2=0, color='grey',
                   edgecolor='k', lw=1)

ax[1].set_ylim(0, emax)
ax[1].set_yticks([])
ax[1].set_xticks([])
ax[1].set_xlabel("DOS", fontsize=18)

plt.savefig('phonon_task4.png')
plt.tight_layout()
plt.savefig('phononSi')
plt.show()