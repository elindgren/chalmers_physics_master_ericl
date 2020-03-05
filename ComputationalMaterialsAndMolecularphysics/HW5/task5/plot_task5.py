# External imports
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
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
# fig, ax = plt.subplots(figsize=(8,6))

fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111, projection='3d')
dbList = list(eosDB.select())

for row in dbList:
    atoms = row.toatoms()
    N = len(atoms.positions)
    dos = row.data['DOS']
    e = row.data['energy'] - row.data['fermi']
    ax.plot(e, N*np.ones(len(e)), dos, alpha=0.8)

ax.set_xlabel(r'Energy (eV)')
ax.set_ylabel(r'Cluster size $N$')
ax.set_zlabel('DOS')
ax.grid()
# plt.tight_layout()
plt.show()