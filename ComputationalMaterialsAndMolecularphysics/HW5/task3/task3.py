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
from ase.vibrations import Vibrations


# Set plot params
plt.rc('font', size=18)          # controls default text sizes
plt.rc('axes', titlesize=18)     # fontsize of the axes title
plt.rc('axes', labelsize=18)     # fontsize of the x and y labels
plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
plt.rc('ytick', labelsize=18)    # fontsize of the tick labels
plt.rc('legend', fontsize=18)    # legend fontsize

def out(s, f):
    ''' Prints s to console and file f. '''
    print(s)
    print(s, file=f)

# Load the database, row by row, and calculate the cohesive energy.
# The cohesive energy is the same as the potential energy since that is
# the energy required to separate the atoms. 
db = connect('../CourseGitRepo/HA5_Al-clusters-initial.db')
vibDB = connect('./vib.db', append=False)
# Sort the clusters based on number of atoms
allClust = list(db.select())
sort = np.argsort([len(clust.numbers) for clust in allClust])
allClust = np.array(allClust)[sort]

for clust in allClust:
    # General info
    atoms = clust.toatoms()
    N = len(atoms.positions)
    # Define calculator - mishin
    mishin = EAM(potential='../CourseGitRepo/HA5_al_potential.alloy') # Set up EAM
    atoms.set_calculator(mishin)
    # Get vibrational spectrum
    v = Vibrations(atoms, name=f'./vibs/vib_{N}')
    str1 = f'--------        N={N}'
    str2 = '        ---------'
    print(str1 +  str2.ljust(40-len(str1)))
    v.run()
    # Get frequencies and DOS - i.e # of states per frequency
    (freq, counts) = np.unique(v.get_frequencies(), return_counts=True)
    freq = np.array(freq)
    dos = np.array(counts)
    # Save to db
    vibDB.write(atoms, data={'frequency': freq, 'DOS': dos})
print('Sucessfully saved all data to DB')


