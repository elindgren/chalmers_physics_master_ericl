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
line_cycler = cycler('linestyle',['-','--',':','-.', '-', '--']) + cycler('color', ['r', 'g', 'b', 'k', 'c', 'y'])
plt.rc('axes', prop_cycle=line_cycler)


N = 6**3 * 4
mAtom = 39.948 # u
aMass = 1.66e-27 # kg
PP = {
    'dt': [],
    'xy': [],
    'xz': [],
    'yz': [],
    'avg': []
}
# Read heatflux
with open('P0Pt.dat', 'r') as f:
    for i, row in enumerate(f):
        if(i>3):
            c = row.rstrip()
            rd = c.split(' ')
            dt = float(rd[1])
            f1 = float(rd[3])
            f2 = float(rd[4])
            f3 = float(rd[5])
            favg = np.average([f1,f2,f3])
            # Save
            PP['dt'].append(dt)
            PP['xy'].append(f1)
            PP['xz'].append(f2)
            PP['yz'].append(f3)
            PP['avg'].append(favg)                                          
for key in PP:
    PP[key] = np.array( PP[key] )

# read volume
V = []
with open('vol.output', 'r') as f:
    # Only final volume is interesting due to running average
    for i, row in enumerate(f):
        if(i>1):
            V.append(float(row.rstrip().split(' ')[1]))
V = np.array( V )

# read pressure
P = []
with open('press.output', 'r') as f:
    # Only final volume is interesting due to running average
    for i, row in enumerate(f):
        if(i>1):
            P.append(float(row.rstrip().split(' ')[1]))
P = np.array( P )


# Plot heatflux correlation functions

# Calculate integral of heatflux function
PP_trap = []
for i in range(len(PP['dt'])):
    dt = PP['dt'][:i]
    PP_avg = PP['avg'][:i]
    PP_trap.append( np.trapz(PP_avg, dt) )

# Plot --------------------
fig, ax1 = plt.subplots(figsize=(8,6))
ax2 = ax1.twinx()
ax1.plot(V, linewidth=2, c='b', linestyle='-', label='V')
ax1.set_ylabel(r'Volume ($\rm Å^3$)', color='b')
ax1.set_xlabel(r'Dump iteration')

ax2.plot(P, linewidth=2, c='r', linestyle='--', label='P')
ax2.set_ylabel(r'Pressure (bar)', color='r')

h1, l1 = ax1.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
handles = [h1[0], h2[0]]
labels = [l1[0], l2[0]]
ax1.legend(handles, labels, loc='upper right')

plt.tight_layout()

#---------------------------

fig, ax1 = plt.subplots(figsize=(8,6))
ax2 = ax1.twinx()
ax1.plot(PP['dt'], PP['avg'], linewidth=2, c='b', linestyle='-', label=r'$\left< P(t)P(0) \right>$')
ax1.set_ylabel(r'$\left< P(t)P(0) \right>$', color='b')
ax1.set_xlabel(r'$t$, ps')

ax2.plot(PP['dt'], PP_trap, linewidth=2, c='r', linestyle='--', label=r'$\int_0^t{\left< P(t)P(0) \right>}$')
ax2.set_ylabel(r'$\int_0^t{\left< P(t)P(0) \right>}$', color='r')

h1, l1 = ax1.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
handles = [h1[0], h2[0]]
labels = [l1[0], l2[0]]
ax1.legend(handles, labels, loc='center right')

plt.tight_layout()

#---------------------------

plt.show()
