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
from ase.phonons import Phonons


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
    str1 = f'--------        N={N}'
    str2 = '        ---------'
    print(str1 +  str2.ljust(40-len(str1)))

    ##### Using Vibrations module
    v = Vibrations(atoms, name=f'./vibs/vib_{N}')
    if not p.isfile(f'./vibs/vib_{N}.all.pckl'):
        print('Running vibration calculation')
        v.run()
        v.combine()  # Combine pickle files
    # Get frequencies and DOS - i.e # of states per frequency
    all_freq = v.get_frequencies()
    # if N==38:
    #     print(v.summary())
    #     print(all_freq)
    (freq, counts) = np.unique(all_freq, return_counts=True)
    fold_freq = v.fold(np.real(freq), np.real(counts), start=0, end=np.real(freq.max()), width=12, normalize=False)
    f_freq = np.array(fold_freq[0])
    f_dos = np.array(fold_freq[1])
    freq = np.array(freq)
    dos = np.array(counts)
    
    # Get number of modes
    # e = False
    # i=0
    # modes = 0
    # while not e:
    #     try:
    #         v.get_mode(i)
    #         modes += 1
    #         i += 1
    #     except:
    #         e = True

    # Save to db
    vibDB.write(atoms, data={'frequency': freq, 'DOS': dos, 'f_freq': f_freq, 'f_dos': f_dos})

    ##### Using Phonons module
    # ph = Phonons(atoms, mishin, delta=0.05, name='./phonons/ph_Si')
    # ph.run()  
    # # ph.combine()
    # ph.read(acoustic=True)
    # Pdos = ph.get_dos(kpts=(20, 20, 20)).sample_grid(npts=60, width=1e-3)
    # print(Pdos)

print('Sucessfully saved all data to DB')


