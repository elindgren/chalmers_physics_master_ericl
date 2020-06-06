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

#! Fixed density and hence V
timestep = 0.01
N = 6**3 * 4
mAtom = 39.948 # u
aMass = 1.66e-27 # kg
kB = 1.38064852e-23 # m2*kg*s-2*K-1
JJ = {
    'dt': [],
    'x': [],
    'y': [],
    'z': [],
    'avg': []
}
# Read heatflux
with open('J0Jt.dat', 'r') as f:
    for i, row in enumerate(f):
        if(i>3):
            c = row.rstrip()
            rd = c.split(' ')
            dt = float(rd[1]) * timestep
            f1 = float(rd[3])
            f2 = float(rd[4])
            f3 = float(rd[5])
            favg = np.average([f1,f2,f3])
            # Save
            JJ['dt'].append(dt)
            JJ['x'].append(f1)
            JJ['y'].append(f2)
            JJ['z'].append(f3)
            JJ['avg'].append(favg)                                          
for key in JJ:
    JJ[key] = np.array( JJ[key] )

# read volume
V = []
with open('vol.output', 'r') as f:
    # Only final volume is interesting due to running average
    for i, row in enumerate(f):
        if(i>1):
            V.append(float(row.rstrip().split(' ')[1]))
V = np.array( V )
rho = N*mAtom*aMass / V[-1] * 1e3/(1e-8)**3 # units: g/cm^3
print(f'Density: rho = {rho:e} g/cm^3') 

# read temperature
T = []
with open('temp.output', 'r') as f:
    # Only final volume is interesting due to running average
    for i, row in enumerate(f):
        if(i>1):
            T.append(float(row.rstrip().split(' ')[1]))
T = np.array( T )
print(f'Temerature: T = {T[-1]:e} K') 


# # read pressure
# P = []
# with open('press.output', 'r') as f:
#     # Only final volume is interesting due to running average
#     for i, row in enumerate(f):
#         if(i>1):
#             P.append(float(row.rstrip().split(' ')[1]))
# P = np.array( P )
# print(f'Pressure: P = {P[-1]:.3f} bar') 


# Calculate integral of heatflux function
JJ_trap = []
for i in range(len(JJ['dt'])):
    dt = JJ['dt'][:i]
    PP_avg = JJ['avg'][:i]
    JJ_trap.append( np.trapz(PP_avg, dt) )
JJ_trap = np.array(JJ_trap) # Å^2*eV^2/ps
T_f = T[-1]
V_f = V[-1] # Å^3
L_qq = JJ_trap * 1/(3*kB*V_f) * (1.602e-19)**2 / 1e-10 # eV^2/(ps*Å*J/K) = (1.609e-19)^2/1e-10 J*K/s*m
ls = L_qq / T_f**2 # W/m*k
# Calculate the thermal conductitivity
l = L_qq[40]
print(f'The thermal conductivity is: {l:e} W/m*K')

# Plot --------------------
fig, ax1 = plt.subplots(figsize=(8,6))
ax2 = ax1.twinx()
ax1.plot(T, linewidth=2, c='b', linestyle='-', label='T')
ax1.set_ylabel(r'Temperature $T$, (K)', color='b')
ax1.set_xlabel(r'Dump iteration')

ax2.plot(V, linewidth=2, c='r', linestyle='--', label='P')
ax2.set_ylabel(r'Volume $V$, (Å^3)', color='r')

h1, l1 = ax1.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
handles = [h1[0], h2[0]]
labels = [l1[0], l2[0]]
ax1.legend(handles, labels, loc='upper right')

plt.tight_layout()
plt.savefig('temp_vol.png')

#---------------------------

fig, ax1 = plt.subplots(figsize=(8,6))
ax2 = ax1.twinx()
ax1.plot(JJ['dt'], JJ['avg'], linewidth=2, c='b', linestyle='-', label=r'$\left< \sigma(t)\sigma(0) \right>$')
ax1.set_ylabel(r'$\left< J(s)J(0) \right>$, ($\rm Å^2eV^2/ps^2$)', color='b')
ax1.set_xlabel(r'$t$, ps')

ax2.plot(JJ['dt'], L_qq, linewidth=2, c='r', linestyle='--', label=r'$\eta$')
ax2.set_ylabel(r'$\lambda = \frac{1}{3kBVT^2}\int_0^t{\left< J(s)J(0) \right>}ds$, (W/m/K)', color='r')

h1, l1 = ax1.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
handles = [h1[0], h2[0]]
labels = [l1[0], l2[0]]
ax1.legend(handles, labels, loc='center right')
ax1.set_xlim(0,2)
plt.tight_layout()
plt.savefig('auto_thremcond.png')

#---------------------------

# plt.show()
