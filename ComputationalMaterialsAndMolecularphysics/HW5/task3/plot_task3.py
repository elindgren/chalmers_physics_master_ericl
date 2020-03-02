# External imports
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

# ASE
from ase import Atoms
from ase.db import connect


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

# Read DB
vibDB = connect('./vib.db')

# Plot vibrational spectra
vFile = open('vib.txt', 'w')
out(s='-----------------    Cluster Vibrations    -----------------', f=vFile)
fig, ax = plt.subplots(figsize=(8,6))
for clust in vibDB.select():
    atoms = clust.toatoms()
    freq = clust.data['frequency']
    dos = clust.data['DOS']
    # Extract info
    N = len(atoms.positions)
    nModes = len(freq)
    isC = np.iscomplex(freq).sum()
    # Print info
    str1 = f'| Al{N} '
    str2 = 'Modes: ' 
    str3 = f' {nModes} '
    str4 = '    Imaginary frequencies: '
    str5 = f' {isC}  |'
    out(str1 + str2.rjust(20-len(str1)) +  str3.rjust(15-len(str2)) + str4 + str5, f=vFile)
    # Plot
    freq = np.real(freq)
    ax.plot(freq, dos, linewidth=2, markersize=2, marker='o', alpha=0.7, label=f'Al{N}')
ax.grid()
ax.legend(loc='best')
ax.set_xlabel(r'Frequency ($\rm cm^{-1}$)')
ax.set_ylabel('DOS')
plt.tight_layout()
plt.savefig('figTask3.png')
out(s='------------------------------------------------------------', f=vFile)