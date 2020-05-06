# External imports
import scipy as sc
import numpy as np
import matplotlib.pyplot as plt

# Set plot params
plt.rc('font', size=18)          # controls default text sizes
plt.rc('axes', titlesize=18)     # fontsize of the axes title
plt.rc('axes', labelsize=18)     # fontsize of the x and y labels
plt.rc('xtick', labelsize=16)    # fontsize of the tick labels
plt.rc('ytick', labelsize=16)    # fontsize of the tick labels
plt.rc('legend', fontsize=14)    # legend fontsize

'''
    This script calculates the diffusion activiation energy in a 4x4x4 supercell of Fe. 
    
    Devloped by Eric Lindgren, SNUID: 2020-81634
    April 2020
'''

dt = 0.002 # timestep, ps
T = []  # temperatures
tmax = []
MSD_sq = []  # absolute mean squared displacements squared
D = []
f = 'output.Fe-npt-msd'

with open(f, 'r') as f:
    f_cont = f.readlines()
    
for i, row in enumerate(f_cont):
    if '----STARTING----' in row:
        i_s = i+7
    elif '----FINISHED' in row:
        d = f_cont[i_s:i+2]
        d.pop(-2)
        # Extract tmax and temperature
        last_it = list(filter( None, d[-1].rstrip().split(' ') ))
        tm = float( last_it[0] )*dt
        tmax.append(tm)
        T.append(float( last_it[6] ))
        # Extract this MSD
        m_sq = [float(list(filter( None, d_i.rstrip().split(' ') ))[-1]) for d_i in d]
        MSD_sq.append(m_sq)
        # Calculate diffusion coefficient
        t = np.linspace(0, tm, len(m_sq))
        D_T = np.polyfit(t, m_sq, deg=1)[0]/6
        D.append(D_T)
    

T = np.array(T)
D = np.array(D)
MSD_sq = np.array(MSD_sq)

# Perform Arrhenius plot and calculate D0 and E
E_R, lnD0 = np.polyfit(1/T, np.log(D), deg=1) 
kB = 8.61733e-5 # eV/K
Na = 6.02214e23 # mol^-1
R = kB*Na # eV/(K*mol)
print(f'D0: {np.exp(lnD0):.3f} ' + r'Å^2/ps, ' + f'{(np.exp(lnD0)*(1e-10)**2/1e-12/1e-12):.3f}e-12 m^2/s')
print(f'E: {-E_R*R/Na:.3f} ' + r'eV per atom, ' + f'{-E_R*R*1.602e-19/1e3:.3f} kJ/mol')
fig, ax = plt.subplots(figsize=(8,6))
ax.scatter(1/T, np.log(D), label='Raw data')
ax.plot(1/T, 1/T*E_R+lnD0, linewidth=2, label='Linear fit')
ax.set_xlabel(r'1/T, K')
ax.set_ylabel(r'$\rm ln\left( D(T) \right)$, $\rm ln \left(  Å^2/ps \right)$')
ax.grid()
ax.legend(loc='best')
plt.tight_layout()
plt.savefig('arrhenius.png')


# Plot absolute MSD
fig, ax = plt.subplots(figsize=(10,9))
for i, temp in enumerate(T):
    msd_sq = MSD_sq[i]
    t = np.linspace(0, tmax[i], len(msd_sq))    
    ax.scatter(t, msd_sq, color=f'C{i}', s=12, alpha=0.7, label=f'Raw data, T={temp:.3f} K')
    ax.plot(t, t*6*D[i], color=f'C{i}', linewidth=2, alpha=1, label=f'Linear fit, T={temp:.3f} K')
ax.set_xlabel(r'$t$, ps')
ax.set_ylabel(r'MSD, $\rm Å^2$')
ax.grid()
ax.legend(loc='best')
plt.tight_layout()
plt.savefig('msd_fit.png')