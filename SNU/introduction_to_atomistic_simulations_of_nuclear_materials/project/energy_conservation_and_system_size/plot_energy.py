# Internal imports
from cycler import cycler
import os

# External imports
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Set plot params
plt.rc('font', size=16)          # controls default text sizes
plt.rc('axes', titlesize=16)     # fontsize of the axes title
plt.rc('axes', labelsize=16)     # fontsize of the x and y labels
plt.rc('xtick', labelsize=16)    # fontsize of the tick labels
plt.rc('ytick', labelsize=16)    # fontsize of the tick labels
plt.rc('legend', fontsize=14)    # legend fontsize
line_cycler = cycler('linestyle',['-','--',':','-.', '-', '--']) + cycler('color', ['b', 'r', 'g', 'k', 'c', 'y'])
plt.rc('axes', prop_cycle=line_cycler)

Natoms = 40000

for e_fold in os.listdir('.'):
    if 'energies' in e_fold:
        E = e_fold.split('_')[1]
        print(f'PKA energy: {E}')
        fig, ax = plt.subplots(figsize=(8,6))
        for efile in os.listdir(e_fold):
            tmp = efile.split('_')[0]
            dt = tmp.split('=')[1]
            df = pd.read_csv(f'./{e_fold}/{efile}', names=['time','ke', 'pe', 'etotal'], skiprows=1 )
            ke = df['ke']
            pe = df['pe']
            t = df['time']
            etotal = df['etotal']
            delta_etotal = etotal-etotal[0]
            avg_deltae = np.average(delta_etotal)
            print(f'\tdt={dt} --- Average energy deviation: {avg_deltae:e} eV')
            ax.plot(t, delta_etotal, label=r'$\Delta t$='+f'{dt} ps')
        ax.grid()
        ax.legend(loc='best')
        ax.set_ylabel(r'$\Delta E_{total}$, (eV)')
        ax.set_xlabel('Time, (ps)')
        ax.set_title(f'PKA Energy: {E}, {Natoms} W atoms')
        plt.tight_layout()
        plt.savefig(f'pka_{E}.png')

# TL;DR: Sweetspot: PKA energy of 2500 eV, with a timestep (dynamic) of minimum 0.0002 ps. Took around 1.5 hrs for 0.3 ps of simulation. 
# However, this is only for a PKA with direction 100 (just travels through)