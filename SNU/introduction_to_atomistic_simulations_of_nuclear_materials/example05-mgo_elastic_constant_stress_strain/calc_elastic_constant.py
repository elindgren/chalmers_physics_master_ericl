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
T = 0
a0 = 4.21079
V = a0**3
d_file = 'd.mgo'
px_file = 'px.mgo'
py_file = 'py.mgo'
pz_file = 'pz.mgo'


# Extract strain and E values
d = extract(d_file)
px = extract(px_file)
py = extract(py_file)
pz = extract(pz_file)

# Extract strain
e = (d-a0)/a0

# Convert to numpy arrays and normalize
px = np.array(px)
py = np.array(py)
pz = np.array(pz)
d = np.array(d)

# Perform linefit
f1 = np.polyfit( d, px, deg=1 ) 
C11 = np.polyfit( d, px, deg=1 )[1]  # bar - strain is dimensionless
C12 = np.polyfit( d, py, deg=1 )[1] #! Shouldn't be 1 here, should be 0 - Something wrong in the definition of np.polyfit
C13 = np.polyfit( d, pz, deg=1 )[1]
# Calculate elastic constant
C11 *= 1e5/1e9 # Convert to GPa
C12 *= 1e5/1e9
C13 *= 1e5/1e9
print(f'Obtained elastic constant at T={T} K: \t C11 = {C11:.3f} GPa, C12 = {C12:.3f} GPa, C13 = {C13:.3f} GPa')


# Plot stress-strain
fig, ax = plt.subplots(figsize=(8,6))
ax.plot( e, px, linestyle='-', label=r'$P_{xx}$' )
ax.plot( e, py, linestyle='--', label=r'$P_{yy}$' )
ax.plot( e, pz, linestyle='-.', label=r'$P_{zz}$' )
ax.set_xlabel(r'Strain $e$')
ax.set_ylabel(r'Energy $E$, (eV)')
ax.grid()
ax.legend(loc='best')
plt.tight_layout()
plt.savefig('stress_strain.png')
plt.show()
