# External imports
import numpy as np
import matplotlib.pyplot as plt


# Set plot params
plt.rc('font', size=18)          # controls default text sizes
plt.rc('axes', titlesize=18)     # fontsize of the axes title
plt.rc('axes', labelsize=18)     # fontsize of the x and y labels
plt.rc('xtick', labelsize=14)    # fontsize of the tick labels
plt.rc('ytick', labelsize=18)    # fontsize of the tick labels
plt.rc('legend', fontsize=18)    # legend fontsize


'''
    This script calculates the bulk modulus K given V and P written 
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
v_file = 'v.fe'
p_file = 'p.fe'

# Extract V and P values
V = extract(v_file)
P = extract(p_file)

# Perform linefit
dPdV = np.polyfit( V, P, deg=1 )[0]  # Bar/Å^3

# Calculate bulk modulus
dPdV *= 1e5/1e9  # Convert to GPa/Å^3
V0 = V[2]  # Å^3
K = -V0 * dPdV
print(f'Obtained bulk modulus at T={T} K: \t K = {K:.3f} GPa')


# Plot stress-strain
fig, ax = plt.subplots(figsize=(8,6))
a = V**(1/3)
a0 = a[2]
e = (a-a0)/a0
ax.plot( e, P*1e-4, label=r'$T=$'+f'{T:.2f} K' )
ax.set_xlabel(r'Strain $e$')
ax.set_ylabel(r'Stress $P$, (GPa)')
ax.grid()
ax.legend(loc='best')
plt.tight_layout()
plt.savefig('stress_strain_0K.png')
plt.show()
