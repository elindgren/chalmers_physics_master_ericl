import numpy as np 
from ase.db import connect 
from ase.calculators.eam import EAM 
from ase.io import write 
from ase import Atoms
from ase.build import bulk 
from ase.dft.dos import DOS
from gpaw import GPAW, FermiDirac, PW
from gpaw import restart 


atoms = bulk('Al', crystalstructure='fcc', a=4.043)

k_values = []
Energies = []

for i in range(1,50): 
    calc = GPAW(mode = PW(300), kpts=(i,i,i))
    atoms.set_calculator(calc) 
    k_values.append(i)
    Energy = atoms.get_potential_energy()
    Energies.append(Energy) 

k_values = np.array(k_values)
Energies = np.array(Energies)

np.savetxt(f'Task6_find_k.txt', np.c_[k_values, Energies])
