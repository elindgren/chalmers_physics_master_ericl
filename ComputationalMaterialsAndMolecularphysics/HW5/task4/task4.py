# Internal imports
import os.path as p

# External imports
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

# ASE
from ase import Atoms
from ase.db import connect
from ase.calculators.eam import EAM
from ase.build import bulk
from ase.vibrations import Vibrations


# Set plot params
plt.rc('font', size=18)          # controls default text sizes
plt.rc('axes', titlesize=18)     # fontsize of the axes title
plt.rc('axes', labelsize=18)     # fontsize of the x and y labels
plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
plt.rc('ytick', labelsize=18)    # fontsize of the tick labels
plt.rc('legend', fontsize=18)    # legend fontsize


# DBs
vibDB = connect('./vib.db', append=False)  # DB for vibration spectrum

# Using optimal lattice parameter from task 2
a = 4.043  # A
atoms = bulk('Al', 'fcc', a)

# Attach EAM calculator to atoms
mishin = EAM(potential='../CourseGitRepo/HA5_al_potential.alloy') # Set up EAM
atoms.set_calculator(mishin)

# Get vibrational spectrum
v = Vibrations(atoms, name='./vib_bulk')
if not p.isfile('vib_bulk.all.pckl'):
    print('Running vibration calculation')
    v.run()
    v.combine()  # Combine pickle files
# Get frequencies and DOS - i.e # of states per frequency
(freq, counts) = np.unique(v.get_frequencies(), return_counts=True)
freq = np.array(freq)
dos = np.array(counts)
v.summary()
# Save to db
vibDB.write(atoms, data={'frequency': freq, 'DOS': dos})




