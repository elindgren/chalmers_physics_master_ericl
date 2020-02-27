 # Imports
import numpy as np
import matplotlib.pyplot as plt

# Set plot params
plt.rc('font', size=18)          # controls default text sizes
plt.rc('axes', titlesize=18)     # fontsize of the axes title
plt.rc('axes', labelsize=18)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
plt.rc('ytick', labelsize=18)    # fontsize of the tick labels
plt.rc('legend', fontsize=18)    # legend fontsize



# Load data
spec_x = np.loadtxt('Na8_spectrum_x.dat', skiprows=6)
spec_y = np.loadtxt('Na8_spectrum_y.dat', skiprows=6)
spec_z = np.loadtxt('Na8_spectrum_z.dat', skiprows=6)
energies = spec_x[:,0]
spec = np.sqrt(spec_x[:,1]**2 + spec_y[:,1]**2 + spec_z[:,1]**2)

# Plot
task1 = np.loadtxt('../task1/spectrum_w0.06.dat', skiprows=4)

fig, ax = plt.subplots(figsize=(8,6))
ax.plot(task1[:,0], task1[:,1] / max(task1[:,1]), color='C2', linestyle='-', linewidth='2', label='GPAW Casida (T.1)') # Normalize both spectra since arbitrary units
ax.plot(energies, spec / max(spec),  color='C0', linestyle='-', linewidth='2', label='Time propagation (T.4)')
ax.set_xlabel('Energy, (eV)')
ax.set_ylabel(r'Photoabsorption spectrum, arb. units')
ax.set_xlim(0,7)
ax.grid()
ax.legend(loc='best')
plt.tight_layout()
plt.savefig('task4_spectra.png')
plt.show()
