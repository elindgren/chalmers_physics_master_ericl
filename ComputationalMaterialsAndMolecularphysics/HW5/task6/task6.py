# Internal imports
import time
import pickle

# External imports
import numpy as np

# ASE
from ase import Atoms
from ase.db import connect
from ase.build import bulk
from ase.parallel import world
from ase.dft.dos import DOS

# GPAW
from gpaw import GPAW, PW, restart

''' 
Perform GPAW DFT calculation for the electronic structure of 
bulk Al. First find converge number of k-points. Then, converge
density self-consistently.
'''

def convergeK(atoms, tol=1e-4, kstart=4):
    # Converge total energy by increasing k-space sampling until total energy changes by
    # <10^-4 eV. 

    # DBs
    convDB = connect('./bulk.db', append=False)  # DB for electronic spectrum

    k = kstart
    Etot_old = 1
    Etot_new = 2
    E = []
    ks = []
    i = 1
    while np.abs(Etot_new - Etot_old) > tol:
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
        ks.append(k)
        k += 4
        i += 1
    # Save calculator state and write to DB
    # if world.rank == 0:
    convDB.write(atoms, data={'energies': E, 'ks': ks})
    calc.write('kConverge.gpw')
    print('Written to DB')
    return k, calc


# Using optimal lattice parameter from task 2
a = 4.043  # A
atoms = bulk('Al', 'fcc', a)

# Find optimal k parameter
k, calc = convergeK(atoms, tol=1e-4, kstart=4)
print(f'Optimal k-parameter: k={k}')

# Perform a ground state energy calculation to get the ground state density
atoms.get_potential_energy()

# Save the calculator
calc.write('Al_calc.gpw')
# if world.rank == 0:
print('Calculator saved')

#### Electronic band structure
# if world.rank == 0:
print('Electronic structure calculation started')
atoms, calc = restart('Al_calc.gpw')
# kpts = {'size': (60,60,60), 'path': 'GXWKGLUWLK,UX'}
kpts = {'path': 'GXWKGLUWLK,UX', 'npoints': 60}
calc.set(kpts = kpts, fixdensity=True, symmetry='off')

# calc = GPAW(
#     'Al_calc.gpw',
#     nbands=16,                              # Include more bands than convergence since metallic
#     fixdensity=True,                        # Fixate the density
#     symmetry='off',                         # Check all points along the path
#     kpts={'path': 'GXWKL', 'npoints': 60},
#     convergence={'bands': 8},
#     txt='Al_calc.txt'
# )
# calc.get_potential_energy()  # Converge the system
# # if world.rank == 0:
# print('Electronic structure converged')

# Fix the potential
calc.get_potential_energy()

# Get band structure and dos
Ebs = atoms.calc.band_structure()  # Get the band structure

# if world.rank == 0:
print('Electronic band structure calculated')
kpts = {'size': (40,40,40)}
calc.set(
    kpts = kpts, 
    fixdensity=True,
    symmetry='off',  
)
# Fix the potential
calc.get_potential_energy()

# e, dos = calc.get_dos(spin=0, npts=1001, width=0.5)  # Get energy and density of states
dos = DOS(calc, npts=2000, width=0.1)
d = dos.get_dos()
e = dos.get_energies() 
f = calc.get_fermi_level() 
print('Electronic DOS computed')
e_f = calc.get_fermi_level()  
Edos = {
    'e': e, 
    'dos': d,
    'fermi': e_f
} 

# Save results
pickle.dump( Ebs, open( "Ebs.p", "wb" ) )  # Save the electronic band structure
pickle.dump( Edos, open( "Edos.p", "wb" ) )  # Save the electronic DOS
# if world.rank == 0:
print('Electronic structure calculation completed')
