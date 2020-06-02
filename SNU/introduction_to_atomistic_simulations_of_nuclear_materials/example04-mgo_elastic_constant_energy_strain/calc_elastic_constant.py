# External imports
import numpy as np
import matplotlib.pyplot as plt


# Set plot params
plt.rc('font', size=14)          # controls default text sizes
plt.rc('axes', titlesize=14)     # fontsize of the axes title
plt.rc('axes', labelsize=14)     # fontsize of the x and y labels
plt.rc('xtick', labelsize=14)    # fontsize of the tick labels
plt.rc('ytick', labelsize=14)    # fontsize of the tick labels
plt.rc('legend', fontsize=14)    # legend fontsize


'''
    This script calculates the elastic constant from stress-energy relation 
    to files from a LAMMPS calculation.

    Devloped by Eric Lindgren, SNUID: 2020-81634
    April 2020
'''


def extract(f):
    ''' Extracts the quantity of interest from the .fe file f and returns it in a Numoy array '''
    a = []
    with open(f, 'r') as file:
        for row in file:
            if '=' in row:
                a.append( float(row.split('=')[1].rstrip()) )
    a = np.array( a )
    return a


# Setup
eV = 1.602177e-19
T = 0
d_file = 'd.mgo'
e_file = 'e.mgo'
a0 = 4.21079
V = a0**3

# Extract strain and E values
d = extract(d_file)
E = extract(e_file)

# Extract strain
e = (d-a0)/a0

# Convert to numpy arrays and normalize
E = np.array(E) - E[2]
d = np.array(d)

# Perform linefit
dE2dd2 = np.polyfit( d, E, deg=2 )[2]  # eV - strain is dimensionless
# Calculate elastic constant
dE2dd2 *= eV*1e-9 # Convert to GJ
V *= 1e-30  # m3
C11 = 2/V * dE2dd2
print(f'Obtained elastic constant at T={T} K: \t C11 = {C11:.3f} GPa')


# Plot stress-strain
fig, ax = plt.subplots(figsize=(8,6))
ax.plot( e, E, label=r'$T=$'+f'{T:.2f} K' )
ax.set_xlabel(r'Strain $e$')
ax.set_ylabel(r'Energy $E$, (eV)')
ax.grid()
ax.legend(loc='best')
plt.tight_layout()
plt.savefig('energy_strain.png')
plt.show()
