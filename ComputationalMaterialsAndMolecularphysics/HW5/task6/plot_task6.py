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
plt.show()

# Extract and plot band structure and DOS
