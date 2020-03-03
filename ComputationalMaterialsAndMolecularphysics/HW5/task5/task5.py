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
            kpts=(8,8,8),
            random=True,
            txt=f'./gpaw-out/EOS_{N}.txt'
        )  # Use the same calculator as in task6
        atoms.set_calculator(calc)
        pot_e = atoms.get_potential_energy()  # Self-constistently optimize the electron density
        if world.rank == 0:
            print(f'Cluster Al{N} finished potential energy: {pot_e:.2f}')
        
        # Get the electronic DOS after choosing a suitable bandpath
        atoms.calc.set(
            nbands=16,                              # Include more bands than convergence since metallic
            fixdensity=True,                        # Fixate the density
            symmetry='off',                         # Check all points along the path
            kpts={'path': 'GXWKL', 'npoints': 100},
            convergence={'bands': 8}
        )
        e, dos = atoms.calc.get_dos(spin=0, npts=100, width=None)
        e_f = calc.get_fermi_level()  
        e -= e_f  # Subtract the fermi level from the energy    
        end = time.time()


        if world.rank == 0:
            print(f'Cluster Al{N} finished ---- Time: {(end-start):.2f} s')
            eosDB.write(atoms, data={'energy': e, 'DOS': dos})
    else:
        if world.rank == 0:
            print(f'Skipping Al{N}')