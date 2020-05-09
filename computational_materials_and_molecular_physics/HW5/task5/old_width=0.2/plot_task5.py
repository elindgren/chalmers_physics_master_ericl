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

fig = plt.figure(figsize=(9,6))
ax = fig.add_subplot(111, projection='3d')
dbList = list(eosDB.select())

for i,row in enumerate(dbList):
    atoms = row.toatoms()
    N = len(atoms.positions)
    dos = row.data['DOS']
    e = row.data['energy'] - row.data['fermi']
    ax.plot(e, N*np.ones(len(e)), dos, alpha=1-i*0.1, linewidth=2, label=f'Al{N}')

ax.legend(loc='upper left')
ax.set_xlabel(r'Energy relative to $\epsilon_F$ (eV)', labelpad=15)
ax.set_ylabel(r'Cluster size $N$', labelpad=15)
ax.set_zlabel(r'DOS ($\rm m^{-3}J^{-1}$)', labelpad=10) # TODO set proper units
ax.grid()
plt.tight_layout()
plt.savefig('dos_task5')
plt.show()