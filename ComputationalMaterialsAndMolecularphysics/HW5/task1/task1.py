# External imports
import numpy as np
import matplotlib.pyplot as plt

# ASE
from ase import Atoms
from ase.db import connect
from ase.calculators.eam import EAM


# Set plot params
plt.rc('font', size=12)          # controls default text sizes
plt.rc('axes', titlesize=12)     # fontsize of the axes title
plt.rc('axes', labelsize=12)     # fontsize of the x and y labels
plt.rc('xtick', labelsize=12)    # fontsize of the tick labels
plt.rc('ytick', labelsize=12)    # fontsize of the tick labels
plt.rc('legend', fontsize=12)    # legend fontsize

def out(s, f):
    ''' Prints s to console and file f. '''
    print(s)
    print(s, file=f)

# Load the database, row by row, and calculate the cohesive energy.
# The cohesive energy is the same as the potential energy since that is
# the energy required to separate the atoms. 
db = connect('../CourseGitRepo/HA5_Al-clusters-initial.db')
# Sort the clusters based on number of atoms
allClust = list(db.select())
sort = np.argsort([len(clust.numbers) for clust in allClust])
allClust = np.array(allClust)[sort]

Efile = open('ECoh.txt', 'w')
out(s='-----------    Cluster energies    ----------', f=Efile)
for clust in allClust:
    # General info
    atoms = clust.toatoms()
    N = len(atoms.positions)
    # Calculate cohesive energy
    mishin = EAM(potential='../CourseGitRepo/HA5_al_potential.alloy') # Set up EAM
    atoms.set_calculator(mishin)
    ECoh = np.abs(atoms.get_potential_energy() / N)  # Cohesive energy is the positive potential energy
    
    # Print results
    str1 = f'| N={N}'
    str2 = 'Cohesive energy/atom:' 
    str3 = f' {ECoh:.4f} eV   |'
    out(str1 + str2.rjust(30-len(str1)) +  str3.rjust(35-len(str2)), f=Efile)
out('---------------------------------------------', f=Efile)

# Plot and save a picture of the potential
mishin.plot()
plt.tight_layout()
plt.savefig('al_potential.png')
# plt.show()

# The potential may not be perfectly suited for the smaller nanoparticles, since they don't exhibit the same symmetries as the bulk Al. 
# Bulk Al is crystalline fcc.