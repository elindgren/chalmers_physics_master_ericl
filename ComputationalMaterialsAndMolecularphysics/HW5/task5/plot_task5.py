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
eosDB = connect('./eos.db')

# Extract and plot convergence data
fig, ax = plt.subplots(figsize=(8,6))
dbList = list(eosDB.select())

dos = dbList[3].data['DOS']
e = dbList[3].data['energy']
ax.plot(e, dos)

ax.set_xlabel(r'Energy (eV)')
ax.set_ylabel('DOS')
ax.grid()
plt.tight_layout()
plt.show()