# ASE
from ase import Atoms
from ase.dft.bandgap import bandgap
from ase.parallel import world

# GPAW
from gpaw import GPAW, restart


# Restart electronicSi calculator and calculate bandgap

_, calc = restart('Si_electrons.gpw')

# Indirect bandgap
gap, p1, p2 = bandgap(calc, direct=False, output='indirectBandgap.txt')
print(f'Indirect bandgap: {gap:.2f} eV')
gap, p1, p2 = bandgap(calc, direct=True, output='directBandgap.txt')
print(f'Direct bandgap: {gap:.2f} eV')