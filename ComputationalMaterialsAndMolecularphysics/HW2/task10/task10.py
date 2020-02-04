# Built-in packages
import time

# Third-party packages
import numpy as np

from ase import Atoms
from ase.ga.utilities import closest_distances_generator
from ase.visualize import view
from ase.optimize import GPMin
from ase.db import connect
from ase.io import write

from gpaw import GPAW, FermiDirac

def get_most_stable_cluster(db):
    '''Fetches the most stable cluster in the db'''
    id_min, E_min = 0, 0
    for row in db.select(relaxed=True):
        # Find lowest energy candidate
        if row.energy < E_min:
            E_min = row.energy
            id_min = row.id
    print(f'Minimum energy for Na{N} cluster: e = {E_min:.4f} eV with id={id_min}')
    a = db.get(f'id={id_min}').toatoms()
    return a

def calculate_stable_wavefunction(N, ef, view=False):
    print(f'-------------- \t Na{N} \t --------------')
    start = time.time()
    #**** Initialize system ****#
    dirpath_t1 = f'../task1/Na{N}/'  # Path to the results from task 1
    db = connect(f'{dirpath_t1}gadb.db')
    stable_clust = get_most_stable_cluster(db)
    #**** Define the calculator ****#
    calc = GPAW(nbands=10,
                h=0.25,
                txt=f'Na{N}_out.txt',
                occupations=FermiDirac(0.05),
                setups={'Na': '1'},
                mode='lcao',
                basis='dzp')
    #**** Relax the system ****#
    stable_clust.set_calculator(calc)
    dyn = GPMin(stable_clust, trajectory=f'Na{N}_relax_clust.traj', logfile=f'Na{N}_relax_clust.log')
    print(f'****\tRelaxing system of {N} atoms\t****')
    dyn.run(fmax=0.02, steps=100)
    #**** Calculate energy and wavefunction ****#
    e = stable_clust.get_potential_energy()  # Note opposite signa from ga.py
    print(f'Na{N} stable energy: {e} eV', file=ef)

    #**** Write wavefunctions to cube files ****
    basename=f'Na{N}'
    nbands = calc.get_number_of_bands()
    for band in range(nbands):
        wf = calc.get_pseudo_wave_function(band=band)
        fname = '{0}_{1}.cube'.format(basename, band)
        print(f'writing wf {band} to file {fname}----'.rjust(12))
        write(fname, stable_clust, data=wf)
    #**** Finish ****
    end = time.time()
    print(f'****    Elapsed time: {(end-start):.2f} s    ' + '****'.rjust(12))
    print('--------------------------------------------\n')
    #*****************************#

Ns = [6,7,8]
ef = open('energies.txt', 'w+')
for N in Ns:
    calculate_stable_wavefunction(N, ef)
ef.close()
