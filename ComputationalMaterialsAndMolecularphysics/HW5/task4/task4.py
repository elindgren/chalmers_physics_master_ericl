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
from ase.phonons import Phonons
from ase.visualize import view


# Set plot params
plt.rc('font', size=16)          # controls default text sizes
plt.rc('axes', titlesize=16)     # fontsize of the axes title
plt.rc('axes', labelsize=16)     # fontsize of the x and y labels
plt.rc('xtick', labelsize=16)    # fontsize of the tick labels
plt.rc('ytick', labelsize=16)    # fontsize of the tick labels
plt.rc('legend', fontsize=16)    # legend fontsize

# DBs
vibDB = connect('./vib.db', append=False)  # DB for vibration spectrum
phonDB = connect('./ph.db', append=False)  # DB for phonon structure and DOS
# Using optimal lattice parameter from task 2
a = 4.043  # A
atoms = bulk('Al', 'fcc', a) 
# view(atoms)  # DEBUG

#### Calculate vibrational spectrum
# Attach EAM calculator to atoms
mishin = EAM(potential='../CourseGitRepo/HA5_al_potential.alloy') # Set up EAM
atoms.set_calculator(mishin)

# Get vibrational spectrum
v = Vibrations(atoms, name='./vibs/vib_bulk')
v.run()
v.summary()
# Get frequencies and DOS - i.e # of states per frequency
(freq, counts) = np.unique(v.get_frequencies(), return_counts=True)
freq = np.array(freq)
dos = np.array(counts)
# Save to db
vibDB.write(atoms, data={'frequency': freq, 'DOS': dos})

#### Calculate band structure 
N=7 # Use a supercell of size 7x7x7 
ph = Phonons(atoms, mishin, name='./phonons/ph_bulk', supercell=(N,N,N), delta=0.05)
ph.run()
# Read the results from the run and obtain the bandpath and DOS
ph.read(acoustic=True)
# ph.clean()

# lat.plot_bz(show=True) # Visualize Brillouin zone
# Al has Space Group 225 [https://materialsproject.org/materials/mp-134/]
# from which Bilbao Cryst gives us 
# Also here is given optimal vectors https://wiki.fysik.dtu.dk/ase/ase/dft/kpoints.html
# Use default path for now
# And here [https://wiki.fysik.dtu.dk/ase/ase/dft/kpoints.html#ase.dft.band_structure.BandStructure]
path = atoms.cell.bandpath(path='GXWKGLUWLK,UX', density=100)

bs = ph.get_band_structure(path)
print(f'Number of phonon modes: {bs.energies.shape[2]}')
dos = ph.get_dos(kpts=(20,20,20)).sample_grid(npts=100, width=1e-3)

# Plot phonon spectrum and DOS
fig, ax = plt.subplots(1,2, figsize=(12,6))
emax = 0.05
bs.plot(ax=ax[0], emin=0.0, emax=emax)

# ax[1] = fig.add_axes([.8, .07, .17, .85])
ax[1].fill_between(dos.weights[0], dos.energy, y2=0, color='grey',
                   edgecolor='k', lw=1)

ax[1].set_ylim(0, emax)
ax[1].set_yticks([])
ax[1].set_xticks([])
ax[1].set_xlabel(r'DOS ($\rm m^{-3}J^{-1}$)') # TODO set proper units

plt.savefig('phonon_task4.png')
plt.tight_layout()
plt.show()





