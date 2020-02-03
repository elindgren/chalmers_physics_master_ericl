# Built-in packages
import time

# Third-party packages
import numpy as np

from ase import Atoms
from ase.ga.utilities import closest_distances_generator
from ase.visualize import view
from ase.optimize import GPMin
from ase.parallel import paropen
from ase.parallel import parprint
from ase.io import read

from gpaw import GPAW, FermiDirac, PW


np.random.seed(1)
d = 20  # cell size
db_file = 'expect.db'  # Filename for each of the databases

def test_params_Na6(m, f, nb):
    print(f'************ Na6 - Mode: {m} - Functional: {f} - nbands={nbands} ************')
    start = time.time()
    structpath=f'../Na-clusters-GA-search/Na6-structures/'  
    #**** Initialize system ****#
    clust0 = read(filename=f'{structpath}christmas-tree.xyz', format='xyz')
    clust1 = read(filename=f'{structpath}half-decahedron.xyz', format='xyz')
    #**** Define the calculator ****#
    if m=='pw':
        calc = GPAW(nbands=nb,
                    h=0.25,
                    txt=f'outfiles/Na6_{m}_{f}_{nb}_out.txt',
                    occupations=FermiDirac(0.05),
                    setups={'Na': '1'},
                    mode=PW(350))
    elif m=='fd':
        calc = GPAW(nbands=nb,
                    h=0.25,
                    txt=f'outfiles/out.txt',
                    occupations=FermiDirac(0.05),
                    setups={'Na': '1'},
                    mode='fd')
    elif m=='lcao':
        calc = GPAW(nbands=nb,
                    h=0.25,
                    txt=f'outfiles/out.txt',
                    occupations=FermiDirac(0.05),
                    setups={'Na': '1'},
                    mode='lcao',
                    basis='dzp')
    #**** Relax the system ****#

    clust0.set_calculator(calc)
    dyn0 = GPMin(clust0, trajectory=f'trajectories/Na6_0_{m}_{f}_{nb}_relax_clust.traj', logfile=f'logfiles/Na6_0_{m}_{f}_{nb}_relax_clust.log')
    clust1.set_calculator(calc)
    dyn1 = GPMin(clust1, trajectory=f'trajectories/Na6_1_{m}_{f}_{nb}_relax_clust.traj', logfile=f'logfiles/Na6_1_{m}_{f}_{nb}_relax_clust.log')
    
    print(f'****    Relaxing 1st system of atoms    ****')
    dyn0.run(fmax=0.02, steps=1)
    print(f'****    Relaxing 2nd system of atoms    ****')
    dyn1.run(fmax=0.02, steps=1)

#     #**** Calculate energy and wavefunction ****#
    e0 = clust0.get_potential_energy()  # Note opposite signa from ga.py
    e1 = clust1.get_potential_energy()
    e0 =1
    e1 =2
    ef = open(f'Na6_e_table.txt', 'a')
    #**** Print to file *****#
    ef.write('-'*78 +'\n')
    ef.write(f'# Mode: {m}' + '|'.rjust(15-len(f'# Mode: {m}')) + '\n')
    subrow1 = f'# XC: {f}' + '|'.rjust(15-len(f'# XC: {f}'))
    subrow2 = f'{e0:.8f} eV \t\t {e1:.8f} eV'.rjust(50-len(subrow1))
    ef.write(subrow1+subrow2 + '\n')
    ef.write(f'# nbands = {nb}' + '|'.rjust(15-len(f'# nbands = {nb}')) + '\n')
    
    end = time.time()
    print(f'****    Elapsed time: {(end-start):.2f} s    ' + '****'.rjust(12))
    print('*'*77 + '\n')
    #*****************************#

# Parameters to perform grid search over
modes = ['pw', 'fd', 'lcao']
functionals=['LSDA', 'PBE']
nbands = [10, 15]

# Prepare file with headers
ef = open(f'Na6_e_table.txt', 'w')
subheader1 = ' Parameters ' + '|'.rjust(15-len(' Parameters '))
subheader2 = '\t\t\t E0 \t\t\t\t  E1'
ef.write(subheader1+subheader2+'\n')
ef.close()
# Iterate over all parameters
print('*'*50 + '\tStarting calculation\t' + '*'*50)
start = time.time()
for m in modes:
    for f in functionals:
        for nb in nbands:
            test_params_Na6(m, f, nb)
end = time.time()
print('#'*25 + (f'\tTotal time: {(end-start):.2f} s\t' + '#'*25).rjust(20))
ef = open(f'Na6_e_table.txt', 'a')
ef.write('-'*78)
ef.close()
