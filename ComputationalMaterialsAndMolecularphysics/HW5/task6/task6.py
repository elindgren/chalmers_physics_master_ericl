# Internal imports
import time

# External imports
import numpy as np

# ASE
from ase import Atoms
from ase.db import connect
from ase.build import bulk
from ase.parallel import world

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
tol = 1e-3
ks = [4*i for i in np.arange(1,11)]  # Nbr of k-points
Etot_old = 1
Etot_new = 2
E = []
i = 1
for k in ks:
    start = time.time()
    Etot_old = Etot_new
    # if world.rank == 0:
    print(f'---- Iteration: {i} ---- k={k} ----')

    calc = GPAW(
            mode=PW(300),                 # cutoff
            kpts=(k, k, k),               # k-points
            txt=f'./gpaw-out/k={k}.txt'   # output file
        )  
    atoms.set_calculator(calc)
    Etot_new = atoms.get_potential_energy()  # Calculates the total DFT energy of the bulk material
    end = time.time()

    # if world.rank == 0:
    print(f'Energy: {Etot_new:.4f} eV ---- Time: {(end-start):.2f} s')
    E.append(Etot_new)
    if np.abs(Etot_new - Etot_old) < tol:
        break
    else:
        # if world.rank == 0:
        i += 1

# Save calculator state and write to DB
# if world.rank == 0:
bulkDB.write(atoms, data={'energies': E, 'ks': ks})
calc.write('bulkConverge.gpw')
print('Written to DB')



# Perform self-consistent density calculation using GPAW


# Save to file
