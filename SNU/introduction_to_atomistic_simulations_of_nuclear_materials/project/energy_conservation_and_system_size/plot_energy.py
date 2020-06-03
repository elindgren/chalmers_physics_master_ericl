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




fig, ax = plt.subplots(figsize=(8,6))

for efile in os.listdir('./energies'):
    tmp = efile.split('_')[0]
    dt = tmp.split('=')[1]
    df = pd.read_csv(f'./energies/{efile}', names=['time','ke', 'pe', 'etotal'], skiprows=1 )
    ke = df['ke']
    pe = df['pe']
    t = df['time']
    etotal = df['etotal']
    delta_etotal = etotal-etotal[0]
    avg_deltae = np.average(delta_etotal)
    print(f'dt={dt} --- Average energy deviation: {avg_deltae:.3f} eV')
    ax.plot(t, delta_etotal, label=r'$\Delta t$='+f'{dt} ps')
ax.grid()
ax.legend(loc='best')
ax.set_ylabel(r'$\Delta E_{total}$, (eV)')
ax.set_xlabel('Time, (ps)')
plt.tight_layout()
plt.show()