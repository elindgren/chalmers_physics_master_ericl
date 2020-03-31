# Internal imports
import pickle

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

def soundV(bs, modes):
    # Get K-points and energies
    k = bs.path.kpts
    e = bs.energies
    N = 30
    l = 2
    v = 0
    print(k.shape)
    print(e.shape)
    hbar = 6.582e-16 # eVs
    fig, ax = plt.subplots(figsize=(8,6))
    for i in range(modes):
        omegar = e[0,:N,i] / hbar
        kr = k[:N, i%3]
        if i==1 or i==2:
            # The only acoustice branches which goes to 0 are these
            # The last mode is overlayed on one of the others
            # if i==1:
            i_off = len([i for i,k in enumerate(kr) if k <= 0.0 and i>0])   # offset to skip zero indexes   
            print(kr[i_off:])
            v = np.abs( (omegar[i_off+l] - omegar[i_off+0]) / (kr[i_off+l] - kr[i_off+0]) ) / 1e10 # Å/s
            print(v)
            ax.plot(kr, omegar*hbar, label=f'{i}')
            ax.plot(kr[i_off:i_off+l], omegar[i_off:i_off+l]*hbar, color='r', label=f'{i}')
    ax.legend(loc='best')
    plt.show()

    return None

# Connect to DB
bulkDB = connect('./bulk.db')

# Extract and plot convergence data
fig, ax = plt.subplots(figsize=(8,6))
dbList = list(bulkDB.select())

ks = dbList[0].data['ks']
E = dbList[0].data['energies']
ax.plot(ks, E)
ax.set_xlabel(r'Number of $k$-points')
ax.set_ylabel('Energy (eV)')
ax.grid()
plt.tight_layout()
plt.savefig('convergenceSi')

################# Extract and plot band electronic band structure and DOS
bs = pickle.load(open( "Ebs.p", "rb" ))
d = pickle.load(open( "Edos.p", "rb" ))

fig, ax = plt.subplots(1,2, figsize=(12,6))

# BS
bs.energies = bs.energies - d['fermi']
emax = 10
emin = -14
bs.plot(filename='', ax=ax[0], show=False, emax=emax, emin=emin)
ax[0].set_ylabel(r'Energy relative to $\epsilon_F$ (eV)')
# lims = (bs.energies.min(), bs.energies.max())
# ax[0].set_ylim(lims)

# DOS
# ax[1].plot(d['e']-d['fermi'], d['dos'])
# ax[1].fill_between( d['dos'], d['e']-d['fermi'], y2=0, color='grey',
#                    edgecolor='k', lw=1)
ax[1].plot(d['dos'], d['e']-d['fermi'])
ax[1].set_xticks([])
# ax[1].set_xlabel(r'DOS ($\rm eV^{-1}$)') # TODO set proper units
ax[1].set_xlabel(r'DOS')
ax[1].set_ylabel(r'Energy relative to $\epsilon_F$ (eV)')
ax[1].set_xticks([])
ax[1].set_ylim((emin,emax))
plt.tight_layout()
plt.savefig('electronicSi')

################# Extract and plot phonon band structure and DOS
bs = pickle.load(open( "Pbs.p", "rb" ))
dos = pickle.load(open( "Pdos.p", "rb" ))

modes = bs.energies.shape[2]
print(f'Number of phonon modes: {modes}')
# Compute sound velocity
# v = soundV(bs, modes)

fig, ax = plt.subplots(1,2, figsize=(12,6))
emax = 0.08
bs.plot(ax=ax[0], emin=0, emax=emax)
ax[0].set_ylabel(r'Energy (eV)')

ax[1].fill_between(dos.weights[0], dos.energy, y2=0, color='grey',
                   edgecolor='k', lw=1)

ax[1].set_ylim(0, emax)
ax[1].set_xticks([])
ax[1].set_ylabel(r'Energy (eV)')
# ax[1].set_xlabel(r'DOS ($\rm eV^{-1}$)') # TODO set proper units
ax[1].set_xlabel(r'DOS')
plt.savefig('phonon_task7.png')
plt.tight_layout()
plt.savefig('phononSi')
plt.show()