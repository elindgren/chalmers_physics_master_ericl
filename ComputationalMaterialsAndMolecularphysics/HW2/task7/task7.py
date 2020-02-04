# General imports
import numpy as np
import matplotlib.pyplot as plt

# ASE
from ase.db import connect
from ase.io import write

# Plot details
plt.rc('font', size=18)          # controls default text sizes
plt.rc('axes', titlesize=18)     # fontsize of the axes title
plt.rc('axes', labelsize=18)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
plt.rc('ytick', labelsize=18)    # fontsize of the tick labels
plt.rc('legend', fontsize=18)    # legend fontsize

def second_most_stable(N):
    dirpath_t1 = f'../task1/Na{N}/'
    db = connect(f'{dirpath_t1}gadb.db')
    
    # Find all energies
    E = np.sort(np.array([row.energy for row in db.select(relaxed=True)]))
    # Find second most stable structure - defined by parameter tol
    tol = 0.04 # A difference of 0.04 eV as compared to lowest state is required
    E0 = E[0]
    idx_E1 = np.where(np.abs(E-E[0])>tol)[0][0]
    E1 = E[idx_E1]
    
    print(f'Na{N}: Energy difference between lowest and second lowest energy state: E1-E0={E1-E0:.4f} eV')
    
    # Plot
    fig, ax = plt.subplots(figsize=(8,6))
    ax.plot(E, label=rf'Sorted $\rm Na_{N}$ cluster energies')
    ax.scatter(0, E0, marker='*', s=100, color='k', label=f'E0={E0:.4f} eV')
    ax.scatter(idx_E1, E1, marker='o', s=100, color='C1', label=f'E1={E1:.4f} eV')
    ax.legend()
    ax.set_xlabel("Index sorted by lowest energy")
    ax.set_ylabel(r"$E$ (eV)")
    ax.grid()
    plt.tight_layout()
    plt.savefig(f'Na{N}_traj.png')
    

Ns = [6,7,8]
for N in Ns:
    second_most_stable(N)