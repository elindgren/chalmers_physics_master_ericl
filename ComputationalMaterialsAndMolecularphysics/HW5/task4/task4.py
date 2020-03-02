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
from ase.visualize import view


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
# view(atoms)


#### Calculate vibrational spectrum
# Attach EAM calculator to atoms
mishin = EAM(potential='../CourseGitRepo/HA5_al_potential.alloy') # Set up EAM
atoms.set_calculator(mishin)

# Get vibrational spectrum
v = Vibrations(atoms, name='./vibs/vib_bulk')
v.run()
# Get frequencies and DOS - i.e # of states per frequency
(freq, counts) = np.unique(v.get_frequencies(), return_counts=True)
freq = np.array(freq)
dos = np.array(counts)
# Save to db
vibDB.write(atoms, data={'frequency': freq, 'DOS': dos})

#### Calculate band structure 
lat = atoms.cell.get_bravais_lattice()
print(lat.description())
# lat.plot_bz(show=True) # Visualize Brillouin zone
# Al has Space Group 225 [https://materialsproject.org/materials/mp-134/]
# from which Bilbao Cryst gives us 
# Use default path for now
path = atoms.cell.bandpath(path='GXWKGLUWLK', density=10)

# Set the path to the calculator
mishin.set(kpts=path, symmetry='off', fixdensity=True)  # Fix the electron density

# Start band structure calculation using EAM (fix electron density)
print(atoms.get_potential_energy())
# bs = mishin.band_structure()
# bs.write('bsTask4.json')


