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


# Load data - skip until frequency data
with open('output.phonon', 'r') as f:
    dmp = f.readlines()
    for i, line in enumerate(dmp):
        if 'k-vector(1/A)	frequency(cm-1)' in line:
            fdata = []
            sub_dmp = dmp[i+1:]
            for l in sub_dmp:
                tmp = [float(i) for i in l.rstrip().split('\t') if not i=='']
                fdata.append(tmp)
fdata = np.array( fdata )

# plot
fig, ax = plt.subplots(figsize=(8,6))
ax.plot(fdata[:,0], np.abs(fdata[:,1]), linewidth=2, label='1D chain 8Ni')
ax.set_xlabel('k, (1/Å)')
ax.set_ylabel('frequency, (cm-1)')
ax.legend(loc='best')
ax.grid()
plt.tight_layout()
plt.show()
