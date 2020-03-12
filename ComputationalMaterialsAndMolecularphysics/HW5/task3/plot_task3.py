# External imports
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
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
out(s='-----------------------------    Cluster Vibrations    -----------------------------', f=vFile)
fig, ax = plt.subplots(figsize=(9,6))
for i, clust in enumerate(vibDB.select()):
    atoms = clust.toatoms()
    freq = clust.data['frequency']
    dos = clust.data['DOS']
    # Extract info
    N = len(atoms.positions)
    nModes = len(freq)
    isC = np.iscomplex(freq).sum()
    isZ = len(freq) - np.count_nonzero(np.real(freq))
    # Print info
    str1 = f'| Al{N} '
    str2 = 'Modes: ' 
    str3 = f' {nModes} '
    str4 = '    Imaginary frequencies: '
    str5 = f' {isC} '
    str6 = '    Zero frequencies: '
    str7 = f' {isZ} |'
    out(str1 + str2.rjust(20-len(str1)) +  str3.rjust(15-len(str2)) + str4 + str5 + str6 + str7, f=vFile)
    # Plot
    f_freq = clust.data['f_freq']
    f_dos = clust.data['f_dos']
    f_freq = np.real(f_freq)
    ax.axhline(0.05*i, color='k', alpha=0.2)
    ax.plot(f_freq, f_dos/nModes + 0.05*i, alpha=0.7, linewidth=2, label=f'Al{N}')
# ax.grid()
ax.legend(loc='upper left')
ax.set_xlabel(r'Wavenumber ($\rm cm^{-1}$)')
# ax.set_ylabel(r'DOS per mode ($\rm cm$)')
ax.set_ylabel(r'DOS per mode')
ax.set_yticks([])
# ax.set_zlabel(r'DOS ($\rm cm$)', labelpad=10) # TODO set proper units
# ax.set_zlabel(r'DOS')
plt.tight_layout()
plt.savefig('figTask3.png')
out(s='------------------------------------------------------------------------------------', f=vFile)