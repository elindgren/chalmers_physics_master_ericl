# Internal imports
import os

# External imports
import numpy as np
import matplotlib.pyplot as plt


# Set plot params
plt.rc('font', size=18)          # controls default text sizes
plt.rc('axes', titlesize=18)     # fontsize of the axes title
plt.rc('axes', labelsize=18)     # fontsize of the x and y labels
plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
plt.rc('ytick', labelsize=18)    # fontsize of the tick labels
plt.rc('legend', fontsize=18)    # legend fontsize


'''
    This script calculates the bulk modulus K given V and P written 
    to files from a LAMMPS calculation for multiple temperatures.

    Devloped by Eric Lindgren, SNUID: 2020-81634
    April 2020
'''

def parse_VP(f):
    ''' Extracts all P and V values from this file '''
    T = []
    V = []
    P = []
    with open(f, 'r') as file:
        f_cont = file.readlines()

    for i, row in enumerate(f_cont):
        if 'FINISHED' in row:
            n_pres = int(row.rstrip().split('n_pres=')[1])
            data = f_cont[i+1].rstrip().split(' ')
            data = list(filter(None, data))
            T.append(float( data[-4] ))
            V.append(float( data[-1] )**3)
            P.append(float( data[-2] ))

    return T, V, P, n_pres


def bulk_modulus(V, P):
    ''' Calculates the bulk modulus given vectors V and P '''
    # Perform linefit
    dPdV = np.polyfit( V, P, deg=1 )[0]  # Bar/Å^3

    # Calculate bulk modulus
    dPdV *= 1e5/1e9  # Convert to GPa/^3
    V0 = V[0]  # Å^3
    K = -V0 * dPdV
    return K

# Setup
d = 'results_high_temp/'

Ts, Vs, Ps, n_pres = parse_VP( d+'output.Fe-high-temp' )
n_temp = int( len(Ts) / n_pres )

Ts = np.array( Ts )
Vs = np.array( Vs )
Ps = np.array( Ps )

Ks = [[0, 177.815]]  # Prepend result for 0K

fig, ax = plt.subplots(figsize=(8,6))
for i in range(n_temp):
    l = i*n_pres
    T = np.mean( Ts[l:l+n_pres] )
    V = Vs[l:l+n_pres]
    P = Ps[l:l+n_pres]
    K = bulk_modulus( V, P )
    Ks.append([T, K])
    print(f'Obtained bulk modulus at T={T:.3f} K: \t K = {K:.3f} GPa')
    a = V**(1/3)
    a0 = a[0]
    e = (a-a0)/a0
    ax.plot( e, P*1e-4, label=r'$T=$'+f'{T:.2f} K' )

# Plot stress-strain -- make sure we are in the elastic linear region
ax.set_xlabel(r'Strain $e$')
ax.set_ylabel(r'Stress $P$, (GPa)')
ax.grid()
ax.legend(loc='best')
plt.tight_layout()
plt.savefig('stress_strain_high_T.png')


# Plot bulk modulus vs temperature
Ks = np.array( Ks )
fig, ax = plt.subplots(figsize=(8,6))
ax.plot( Ks[:,0], Ks[:,1], markersize=12, marker='o', linewidth=2 )
ax.grid()
ax.set_xlabel(r'Temperature $T$, (K)')
ax.set_ylabel(r'Bulk modulus $K$, (GPa)')
plt.tight_layout()
plt.savefig('bulk_modulus_high_T.png')

plt.show()