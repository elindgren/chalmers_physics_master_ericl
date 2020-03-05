# Internal imports
import time
import pickle

# External imports
import numpy as np

# ASE
from ase import Atoms
from ase.db import connect
from ase.build import bulk
from ase.phonons import Phonons
from ase.parallel import world

# GPAW
from gpaw import GPAW, PW, restart

'''
Perform DFT calculation for electronic and phononic band structures and density of states for Si.
Uses a GPAW calculator with a PW basis set. Inspiration taken from this example:
https://wiki.fysik.dtu.dk/gpaw/tutorials/bandstructures/bandstructures.html.
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
                mode=PW(200),                 # cutoff - lower for computational efficiency
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

# Define the Si bulk-structure
atoms = bulk('Si', 'diamond', 5.43)
if world.rank == 0:
    print('System created')

# Find optimal k parameter
k, calc = convergeK(atoms, tol=1e-4, kstart=4)
if world.rank == 0:
    print(f'Optimal k-parameter: k={k}')

# Perform a ground state energy calculation to get the ground state density
atoms.get_potential_energy()

# Save the calculator
calc.write('Si_calc.gpw')
if world.rank == 0:
    print('Calculator saved')

#### Electronic band structure
# if world.rank == 0:
# print('Electronic structure calculation started')
# calc = GPAW(
#     'Si_calc.gpw',
#     nbands=16,                              # Include more bands than convergence since metallic
#     fixdensity=True,                        # Fixate the density
#     symmetry='off',                         # Check all points along the path
#     kpts={'path': 'GXWKL', 'npoints': 60},
#     convergence={'bands': 8},
#     txt='Si_calc.txt'
# )
# calc.get_potential_energy()  # Converge the system
# # if world.rank == 0:
# print('Electronic structure converged')

atoms, calc = restart('Si_calc.gpw')
# kpts = {'size': (20,20,20)}
kpts = {'path': 'GXWKL', 'npoints': 60}
calc.set(
    kpts = kpts, 
    fixdensity=True,
    symmetry='off',  
)

# Fix the potential
calc.get_potential_energy()

# Get band structure and dos
Ebs = atoms.calc.band_structure()  # Get the band structure
if world.rank == 0:
    print('Electronic band structure calculated')

# Set new k-mesh to the calculator to get a nice DOS
kpts = {'size': (28,28,28)}
calc.set(
    kpts = kpts, 
    fixdensity=True,
    symmetry='off',  
)
# Fix the potential
calc.get_potential_energy()

e, dos = calc.get_dos(spin=0, npts=1001, width=0.2)  # Get energy and density of states
e_f = calc.get_fermi_level()  
Edos = {
    'e': e, 
    'dos': dos,
    'fermi': e_f
} 

# Save results
pickle.dump( Ebs, open( "Ebs.p", "wb" ) )  # Save the electronic band structure
pickle.dump( Edos, open( "Edos.p", "wb" ) )  # Save the electronic DOS
calc.write('Si_electrons.gpw')
if world.rank == 0:
    print('Electronic structure calculation completed')


#### Phononic band structure
# if world.rank == 0:
print('Phononic structure calculation started')
atoms, calc = restart('Si_calc.gpw')
# kpts = {'size': (20,20,20)}
calc.set(
    symmetry='off',  
)


# Set up the ASE phonon calculator
N = 2  # Use a 2x2x2 supercell
ph = Phonons(atoms, calc, supercell=(N, N, N), delta=0.05, name='./phonons/ph_Si')

# Run the phonon calculation
if world.rank == 0:
    print('******** Phonon calculation started *********')
ph.run() 
if world.rank == 0:
    print('******** Phonon calculation completed *********')
ph.read(acoustic=True)

# Define BZ-path - use the same as for the electronic calculation
path = atoms.cell.bandpath('GXWKL', npoints=60)

# Fetch band structure and dos
if world.rank == 0:
    print('******** Calculating phononic band structure *********')
Pbs = ph.get_band_structure(path)
if world.rank == 0:
    print('******** Phononic band structure calculated *********')
Pdos = ph.get_dos(kpts=(20, 20, 20)).sample_grid(npts=60, width=1e-3)
if world.rank == 0:
    print('******** Phononic DOS calculated *********')

# Save results
pickle.dump( Pbs, open( "Pbs.p", "wb" ) )  # Save the phononic band structure
pickle.dump( Pdos, open( "Pdos.p", "wb" ) )  # Save the phononic DOS
# calc.write('Si_phonons.gpw') # Don't need to save this calc
if world.rank == 0:
    print('Phononic structure calculation completed')