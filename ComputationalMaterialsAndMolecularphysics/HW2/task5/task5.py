# General imports
import numpy as np
import matplotlib.pyplot as plt

# ASE
from ase.db import connect
from ase.io import write
from ase.visualize import view

# Plot details
plt.rc('font', size=18)          # controls default text sizes
plt.rc('axes', titlesize=18)     # fontsize of the axes title
plt.rc('axes', labelsize=18)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
plt.rc('ytick', labelsize=18)    # fontsize of the tick labels
plt.rc('legend', fontsize=18)    # legend fontsize


def get_most_stable_cluster(db):
    id_min, E_min = 0, 0
    for row in db.select(relaxed=True):
        # Find lowest energy candidate
        if row.energy < E_min:
            E_min = row.energy
            id_min = row.id
    print(f'Minimum energy for Na{N} cluster: e = {E_min:.4f} eV with id={id_min}')
    a = db.get(f'id={id_min}').toatoms()
    return a

def view_results_and_to_xyz(N):
    dirpath_t1 = f'../task1/Na{N}/'
    dirpath_t5 = f'.'
    db = connect(f'{dirpath_t1}gadb.db')
    a = get_most_stable_cluster(db)
    #**** View atoms - save to a png file ****#
    view(a)
    write(filename=f'{dirpath_t5}Na{N}_min.png', images=a)
    #**** Save atoms to XYZ ****#
    write(filename=f'{dirpath_t5}Na{N}_min.xyz', images=a)

Ns = [6,7,8]
for N in Ns:
    view_results_and_to_xyz(N)