# Built-in packages
import time

# ASE
from ase import Atoms
from ase.io import read, write
from ase.parallel import world

# GPAW
from gpaw import GPAW


# Load Na8
atoms = read('../Na-tddft/Na8.xyz')
atoms.center(vacuum=8.0) # Add 8 Å of vacuum around the cluster

# Define calculator
calc = GPAW(
    mode = 'fd',
    xc = 'LDA',
    setups = {'Na': '1'},
    h = 0.3,
    nbands = 0,
    txt = 'conv.gpaw-out'
)
atoms.set_calculator(calc)

#------------ Converge ground state ------------
if world.rank==0:
    print(f'------------   Convergence calculation for ground state   ------------')
start = time.time()

atoms.get_potential_energy()
# Save result for ground state
calc.write('groundCalc.gpw', 'all')

#------------ Converge empty states ------------
if world.rank==0:
    print(f'------------   Convergence calculation for empty states   ------------')

calc.set(
    nbands = 110,
    convergence = {'bands': '-10'},
    fixdensity = True,
    eigensolver = 'cg'
)
atoms.get_potential_energy()
# Save results for empty states
calc.write('emptyCalc.gpw', 'all')


end = time.time()
if world.rank==0:
    print('-------- Convergence calculation finished in: ' + f'{(end-start):.2f} s --------'.rjust(34))
    print('----------------------------------------------------------------------')