# External imports
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

# ASE
from ase import Atoms
from ase.db import connect
from ase.calculators.eam import EAM
from ase.build import bulk


# Set plot params
plt.rc('font', size=18)          # controls default text sizes
plt.rc('axes', titlesize=18)     # fontsize of the axes title
plt.rc('axes', labelsize=18)     # fontsize of the x and y labels
plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
plt.rc('ytick', labelsize=18)    # fontsize of the tick labels
plt.rc('legend', fontsize=18)    # legend fontsize


# The true lattice parameter is 4.04 Å - search around this minimum
latParams = np.linspace(3.5, 4.5, 1000)

E = []
for a in tqdm(latParams):
    atoms = bulk('Al', 'fcc', a)

    # Calculate potential energy
    mishin = EAM(potential='../CourseGitRepo/HA5_al_potential.alloy') # Set up EAM
    atoms.set_calculator(mishin)
    E.append(atoms.get_potential_energy())

E = np.abs(E)  # Cohesive energy is abs of potential energy
# Find maxmum cohesive energy and best lattice parameter
maxIdx = np.argmax(E)
aMax = latParams[maxIdx]
EMax = E[maxIdx]
print(f'Optimal lattice parameter: {aMax:.3f} Å, Cohesive energy: {EMax:.3f} eV')

# Plot results
fig, ax = plt.subplots(figsize=(8,6))
ax.plot(latParams, E, linewidth=2, linestyle='-', marker='o', markersize=0, color='C0', label=r'$E_{Coh}$')
ax.scatter(aMax, EMax, marker='*', s=150, c='k',zorder=3, label=r'$a_{optimal}$ ' + f'= {aMax:.3f} Å')
ax.set_xlabel(r'Lattice parameter $a$ (Å)')
ax.set_ylabel(r'Cohesive energy (eV)')
ax.legend(loc='best')
ax.grid()
plt.tight_layout()
plt.savefig('E_task2.png')
# plt.show()

# The lattaice parameter agrees with experimental values up to 2 decimal
# places, and the cohesive energy seems to match those (per atom) 
# from task1 somewhat well. 
# However, the cohesive energy for bulk aluminium averaged over all atoms
# is lower, which could be due to the increased stability of the bulk crystal.
# In the nanoparticles there are atoms on the surface which are bound less
# tight which reduces the average cohesive energy. The effects of these
# ''surface atoms'' decreases with the size of the nanoparticle which can be 
# seen by the cohesive energy increasing with the number of atoms in the 
# nanoparticle. The bulk Al corresponds to the cohesive energy of an infinite
# size nanoparticle, i.e. no edge effects. 

# We expect EAM to give good results, since it is particularly good for FCC metals
# [https://wiki.fysik.dtu.dk/ase/ase/calculators/eam.html].