# Built-in packages
import time

# Third-party packages
import numpy as np

from ase import Atoms
from ase.ga.utilities import closest_distances_generator
from ase.visualize import view
from ase.optimize import GPMin

from gpaw import GPAW, FermiDirac


np.random.seed(1)
d = 10  # cell size
db_file = 'expect.db'  # Filename for each of the databases

def create_relaxed_Na_cluster(N, view=False):
    print(f'************ Na{N} ************')
    start = time.time()
    dirpath=f'./Na{N}/Task3/'  # Path where everything will be saved
    #**** Initialize system ****#
    if N==6:
        clust = Atoms('Na'*6, positions=[(1,1,0),(1,-1,0),(-1,-1,0),(-1,1,0),(0,0,1),(0,0,-1)], cell=(d, d, d))
    else:
        clust = Atoms('Na'*N, positions=[np.random.randn(3) for i in range(N)], cell=(d, d, d))  # random initialization
    clust.center()
    if view:
        view(Na6)
    #**** Define the calculator ****#
    calc = GPAW(nbands=10,
                h=0.25,
                txt=f'{dirpath}out.txt',
                occupations=FermiDirac(0.05),
                setups={'Na': '1'},
                mode='lcao',
                basis='dzp')
    #**** Relax the system ****#
    clust.set_calculator(calc)
    dyn = GPMin(clust, trajectory=f'{dirpath}relax_clust.traj', logfile='{dirpath}relax_clust.log')
    print(f'Relaxing system of {N} atoms')

    end = time.time()
    print(f'****Elapsed time: {end-start} s****')
    print('*****************************\n')
    #*****************************#

Ns = [6,7,8]
for N in Ns:
    create_relaxed_Na_cluster(N)