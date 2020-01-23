import ase
from ase import Atoms
from ase.optimize import BFGS
from ase.calculators.nwchem import NWChem
from ase.io import write
h2 = Atoms('H2', positions=[[0, 0, 0], [0, 0, 0.7]])
h2.calc = NWChem(xc='PBE')
opt = BFGS(h2)
opt.run(fmax=0.02)

write('H2.xyz', h2)
h2.get_potential_energy()