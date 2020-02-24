# Built-in packages
import time

# ASE
from ase import Atoms
from ase.io import read, write
from ase.parallel import world

# GPAW
from gpaw import GPAW
from gpaw.tddft import *

# Load the ground state calculator from task1
td_calc = TDDFT('../task1/groundCalc.gpw')

# Set up time propagation DFT
time_step = 30  # 30 attoseconds
total_time = 45000
iterations = total_time/time_step
kick_strength = 1e-5  # Kick with a field of this strength in x, y and z directions

# Run calculation in each cartesian direction
coord = ['x', 'y', 'z']

if world.rank==0:
    print(f'------------   Task 4 TDDFT started   ------------')
start = time.time()

for i, c in enumerate(coord): 
    if world.rank==0:
        print(f'------------   Propagating {c}-direction   ------------')
    start_c = time.time()
    # Prepare kick
    kick = [0.0]*3
    kick[i] = kick_strength
    # Kick with a delta kick in direction c
    td_calc.absorption_kick(kick_strength=kick)
    # Propagate system and save time dependent dipole moment to 'Na8_dm_{c}.dat' and
    # use 'Na8_td_{c}.dat' as restart file
    td_calc.Propagate(time_step, iterations, f'Na8_dm_{c}.dat', f'Na8_td_{c}.dat')

    # Calculate photoabsorption spectrum and save it 
    photoabsorption_spectrum(f'Na8_dm_{c}.dat', f'Na8_spectrum_{c}.dat')
    end = time.time()

    if world.rank==0:
        print(f'-------- Propagating {c}-direction finished in: ' + f'{(end-start):.2f} s --------'.rjust(34))
    
    
end = time.time()
if world.rank==0:
    print('-------- TDDFT calculation finished in: ' + f'{(end-start):.2f} s --------'.rjust(34))
    print('----------------------------------------------------------------------')
# Load all spectra, take the average and plot in a separate script. 

