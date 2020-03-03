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

# Define the Si bulk-structure
atoms = bulk('Si', 'diamond', 5.43)
if world.rank == 0:
    print('System created')

# Define GPAW calculator - will be used both for electronic and vibrational calculation
calc = GPAW(
    mode=PW(200),
    kpts=(8,8,8),
    random=True,      # Needed to get many electronic bands for our Si
    txt='Si_calc.txt'
)
atoms.set_calculator(calc)

# Perform a ground state energy calculation to get the ground state density
atoms.get_potential_energy()

# Save the calculator
calc.write('Si_calc.gpw')
# if world.rank == 0:
print('Calculator saved')

#### Electronic band structure
# if world.rank == 0:
print('Electronic structure calculation started')
calc = GPAW(
    'Si_calc.gpw',
    nbands=16,                              # Include more bands than convergence since metallic
    fixdensity=True,                        # Fixate the density
    symmetry='off',                         # Check all points along the path
    kpts={'path': 'GXWKL', 'npoints': 60},
    convergence={'bands': 8},
    txt='Si_calc.txt'
)
calc.get_potential_energy()  # Converge the system
# if world.rank == 0:
print('Electronic structure converged')

# Get band structure and dos
Ebs = calc.band_structure()  # Get the band structure
# if world.rank == 0:
print('Electronic band structure calculated')

e, dos = calc.get_dos(spin=0, npts=60, width=0.2)  # Get energy and density of states
e_f = calc.get_fermi_level()  
e -= e_f  # Subtract the fermi level from the energy
Edos = {
    'e': e, 
    'dos': dos
} 

# Save results
pickle.dump( Ebs, open( "Ebs.p", "wb" ) )  # Save the electronic band structure
pickle.dump( Edos, open( "Edos.p", "wb" ) )  # Save the electronic DOS
# if world.rank == 0:
print('Electronic structure calculation completed')

#### Phononic band structure
# if world.rank == 0:
print('Phononic structure calculation started')
calc = GPAW('Si_calc.gpw')  # Load the calculator

# Set up the ASE phonon calculator
N = 7  # Use a 7x7x7 supercell
ph = Phonons(atoms, calc, supercell=(N, N, N), delta=0.05, name='./phonons/ph_Si')

# Run the phonon calculation
ph.run()  
ph.read(acoustic=True)

# Define BZ-path - use the same as for the electronic calculation
path = atoms.cell.bandpath('GXWKL', npoints=60)

# Fetch band structure and dos
Pbs = ph.get_band_structure(path)
Pdos = ph.get_dos(kpts=(20, 20, 20)).sample_grid(npts=100, width=1e-3)

# Save results
pickle.dump( Pbs, open( "Pbs.p", "wb" ) )  # Save the phononic band structure
pickle.dump( Pdos, open( "Pdos.p", "wb" ) )  # Save the phononic DOS
# if world.rank == 0:
print('Phononic structure calculation completed')