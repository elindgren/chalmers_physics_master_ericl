# Built-in packages
import time

# Third-party packages
import numpy as np

from ase import Atoms
from ase.ga.utilities import closest_distances_generator
from ase.visualize import view
from ase.optimize import GPMin

from gpaw import GPAW, FermiDirac

def create_relaxed_Na_cluster(N, view=False):
    print(f'************ Na{N} ************')
    start = time.time()
    #**** Initialize system ****#
    if N==6:
        clust = Atoms('Na'*6, positions=[(1,1,0),(1,-1,0),(-1,-1,0),(-1,1,0),(0,0,1),(0,0,-1)], cell=(d, d, d))
    else:
        clust = Atoms('Na'*N, positions=[np.random.randn(3) for i in range(N)], cell=(d, d, d))  # random initialization
    clust.center()
    if view:
        view(clust)
    #**** Define the calculator ****#
    calc = GPAW(nbands=10,
                h=0.25,
                txt=f'Na{N}_out.txt',
                occupations=FermiDirac(0.05),
                setups={'Na': '1'},
                mode='lcao',
                basis='dzp')
    #**** Relax the system ****#
    clust.set_calculator(calc)
    dyn = GPMin(clust, trajectory=f'Na{N}_relax_clust.traj', logfile='Na{N}_relax_clust.log')
    print(f'**** Relaxing system of {N} atoms ****')
    dyn.run(fmax=0.02, steps=100)
    #**** Calculate energy and wavefunction ****#
    e = clust.get_potential_energy()  # Note opposite signa from ga.py
    e_file = open(f'Na{N}_e_cluster.txt', 'w')
    print(f'Na{N} cluster energy: {e} eV', file=e_file)
    calc.write(f'Na{N}_cluster.gpw', mode='all')
    end = time.time()
    print(f'**** Elapsed time: {end-start} s ****')
    print('*****************************\n')
    #*****************************#

Ns = [6,7,8]
for N in Ns:
    create_relaxed_Na_cluster(N)
