# Internal imports
from cycler import cycler

# External imports
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


'''
    Plot the Buckingham potential model for MgO.

    Developed by Eric Lindgren, SNU ID ericlin
    June 2020
'''

def buckingham(mp, q1, q2, r):
    ''' mp are the model parameters, q1 and q2 are the two charges and r is the interatomic distance '''
    return 14.4*q1*q2/r + mp[0]*np.exp(-r/mp[1]) - mp[2]/r**6


rc = 10.0 # Å, rcut
q = {
    'Mg': 2,
    'O': -2
} # Charge

model = {
    'Mg-Mg' : { 
        'mp' : [0.0, 0.1, 0.0, rc],
        'q1' : q['Mg'],
        'q2' : q['Mg']
    },
    'Mg-O' : {
        'mp' : [821.6, 0.3242, 0.0, rc],
        'q1' : q['Mg'],
        'q2' : q['O']
    },
    'O-O' : {
        'mp' : [22764.0, 0.1490, 27.88, rc],
        'q1' : q['O'],
        'q2' : q['O']
    }
} # Order: A (eV), rho (Å), C (eV Å^6), r_cut (Å)

r = np.linspace(0.1, rc, 100)

fig, ax = plt.subplots(figsize=(8,6))
for l, c in model.items():
    U = buckingham(c['mp'], c['q1'], c['q2'], r)
    ax.plot(r, U, linewidth=2, label=l)
ax.set_xlabel(r'$r_{12}$, (Å)')
ax.set_ylabel('Potential energy, (eV)')
ax.legend(loc='best')
ax.set_ylim(-100,200)
ax.grid()
plt.show()