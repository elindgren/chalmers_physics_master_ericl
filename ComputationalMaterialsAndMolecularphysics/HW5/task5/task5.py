# Internal imports
import time

# External imports
import numpy as np

# ASE
from ase import Atoms
from ase.db import connect
from ase.dft.dos import DOS
from ase.parallel import world

# GPAW
from gpaw import GPAW, PW

''' 
Perform GPAW DFT calculation of the electron density of states
for the nanoparticles with N < 100. 
'''
# Connect to DB
structDB = connect('../CourseGitRepo/HA5_Al-clusters-initial.db')
eosDB = connect('./eos.db', append=False)

# Sort the clusters based on number of atoms
allClust = list(structDB.select())
sort = np.argsort([len(clust.numbers) for clust in allClust])
allClust = np.array(allClust)[sort]

for clust in allClust:
    # General info
    atoms = clust.toatoms()
    N = len(atoms.positions)
    if(N<100):
        start = time.time()
        if world.rank == 0:
            print(f'Calculating EOS for Al{N}')
        # Define electron calculator (GPAW)
        calc = GPAW(
            mode=PW(300),
            txt=f'./gpaw-out/EOS_{N}.txt'
        )  # Use the same calculator as in task6
        atoms.set_calculator(calc)
        pot_e = atoms.get_potential_energy()  # Just do this to connect calculator to cluster
        if world.rank == 0:
            print(f'Cluster Al{N} finished potential energy: {pot_e:.2f}')
        # Calculate DOS using ASE
        # dos = DOS(calc, width=0.2)
        # d = dos.get_dos()
        # e = dos.get_energies()
        e, dos = calc.get_dos(spin=0, npts=201, width=None)
        end = time.time()
        if world.rank == 0:
            print(f'Cluster Al{N} finished ---- Time: {(end-start):.2f} s')
            eosDB.write(atoms, data={'energy': e, 'DOS': dos})
    else:
        if world.rank == 0:
            print(f'Skipping Al{N}')