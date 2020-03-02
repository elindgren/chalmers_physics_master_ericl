# External imports
import numpy as np

# ASE
from ase import Atoms
from ase.db import connect
from ase.dft.dos import DOS

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
        print(f'Calculating EOS for Al{N}')
        # Define electron calculator (GPAW)
        calc = GPAW(
            mode=PW(300),
            xc='PBE',
            out='EOS.txt'
        )  # Use the same calculator as in task6

        # Calculate DOS using ASE
        dos = DOS(calc, width=0.2)
        d = dos.get_dos()
        e = dos.get_energies()
        eosDB.write(atoms, data={'energy': e, 'DOS': dos})
    else:
        print(f'Skipping Al{N}')