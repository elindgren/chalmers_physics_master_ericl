# Internal imports
import os.path as p

# External imports
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

# ASE
from ase import Atoms
from ase.db import connect

# GPAW
from gpaw import GPAW, PW

''' 
Perform GPAW DFT calculation for the electronic structure of 
bulk Al. First find converge number of k-points. Then, converge
density self-consistently.
'''

# DBs
bulkDB = connect('./bulk.db', append=False)  # DB for vibration spectrum

# Using optimal lattice parameter from task 2
a = 4.043  # A
atoms = bulk('Al', 'fcc', a)

# Converge total energy by increasing k-space sampling until total energy changes by
# <10^-4 eV. 
tol = 1e-4
k = 4  # Nbr of k-points
Etot_old = 1
Etot_new = 2

while np.abs(Etot_old - Etot_new ) > tol:
    k *= 1.1  # Increase k by 10%
    calc = GPAW(mode=PW(300),             # cutoff
            kpts=(k, k, k),               # k-points
            txt=f'./gpaw-out/k={k}.txt')  # output file
    atoms.set_calculator(calc)
    

# Perform self-consistent density calculation using GPAW

# Save to file
