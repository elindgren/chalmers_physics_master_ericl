import numpy as np 
from ase.db import connect 
from ase.calculators.eam import EAM 
from ase.io import write 
from ase import Atoms
from ase.build import bulk 
from ase.dft.dos import DOS
from gpaw import GPAW, FermiDirac, PW
from gpaw import restart 
from ase.dft.band_structure import BandStructure

atoms = bulk('Al', crystalstructure='fcc', a=4.043)

k = 46 

calc = GPAW(mode = PW(300))
atoms.set_calculator(calc) 
atoms.get_potential_energy() 
calc.write('Task6.gpw', 'all')

atoms, calc = restart('Task6.gpw')
kpts={'size': (k,k,k)}
calc.set(kpts=kpts, fixdensity=True)
calc.get_potential_energy() 

dos = DOS(calc, npts=2000, width=0.1)
d = dos.get_dos()
e = dos.get_energies() 
f = calc.get_fermi_level() 
d = np.array(d)
e = np.array(e)
f = np.array(f) 

np.savetxt('Task6_finished.txt', np.c_[e, d]) 

calc.set(kpts={'path': 'GXWKGLUWLK,UX', 'npoints': 100}, symmetry='off')
calc.get_potential_energy()
bs = calc.band_structure()
bs.plot(filename='bandstructure.png', show=True, emax=30.0)



