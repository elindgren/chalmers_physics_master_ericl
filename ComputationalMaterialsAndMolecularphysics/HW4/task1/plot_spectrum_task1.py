# Imports
import numpy as np
import matplotlib.pyplot as plt

# Set plot params
plt.rc('font', size=18)          # controls default text sizes
plt.rc('axes', titlesize=18)     # fontsize of the axes title
plt.rc('axes', labelsize=18)     # fontsize of the x and y labels
plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
plt.rc('ytick', labelsize=18)    # fontsize of the tick labels
plt.rc('legend', fontsize=18)    # legend fontsize



# Load data
spec_data = np.loadtxt('spectrum_w0.06.dat', skiprows=4)

# Plot
fig, ax = plt.subplots(figsize=(8,6))
ax.plot(spec_data[:,0], spec_data[:,1], linestyle='-')
ax.set_xlabel('Energy, (eV)')
ax.set_ylabel(r'Photoabsorption spectrum, arb. units')
ax.set_xlim(0,7)
ax.grid()
plt.tight_layout()
plt.savefig('task1_spectra.png')
plt.show()
