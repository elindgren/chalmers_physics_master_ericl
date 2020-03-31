# External imports
import numpy as np
import matplotlib.pyplot as plt
import helper as h


# Set plot params
plt.rc('font', size=18)          # controls default text sizes
plt.rc('axes', titlesize=18)     # fontsize of the axes title
plt.rc('axes', labelsize=18)     # fontsize of the x and y labels
plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
plt.rc('ytick', labelsize=18)    # fontsize of the tick labels
plt.rc('legend', fontsize=18)    # legend fontsize

# Load dumped data from task 1
d = np.load('dumpTask1.npz')
K = d['K_pp'] # K-matrix
omega = d['ediff_p'] # KS eigenvalue difference
n = d['fdiff_p']  # Occupation difference
mux_p = d['mux_p']
muy_p = d['muy_p']
muz_p = d['muz_p']
mu_p = [mux_p, muy_p, muz_p]

# Construct the Omega matrix
Omega = np.diag(v=omega**2, k=0) # First add diagonal
for p in range(len(Omega)):
    for q in range(len(Omega)):
        Omega[p,q] += 2*np.sqrt(n[p]*omega[p]) * K[p,q] * np.sqrt(n[q]*omega[q])

# Obtain egivenalues and eigenvectors from the Omega matrix
eigVal, F = np.linalg.eig(Omega)  # Eigenvalue, eigenvectors
sort = np.argsort(eigVal)
eigVal = eigVal[sort]
F = F[:,sort]

# Extract the excitations in Ha, and convert them to eV
eps = np.sqrt(eigVal) * 27.2

# Calculate the oscillator strength
f = np.zeros(len(eigVal))
for i, e in enumerate(eigVal):
    for alpha, mua_p in enumerate(mu_p):
        f_ia = 2 * np.abs( np.sum( mua_p * np.sqrt(n*omega) * F[:,i] ) )**2
        f[i] += f_ia
    f[i] /= 3 # Average over all dipole moments

# Compare with the discrete spectrum from GPAW
# Plot
disc_data = np.loadtxt('GPAW_discrete.dat', skiprows=0)
fig, ax = plt.subplots(figsize=(8,6))
for i,e in enumerate(eps):
    if i==0:
        # GPAW 
        d_e = disc_data[i,1]
        I_e = disc_data[i,2]
        ax.plot([d_e, d_e], [0, I_e] / max(disc_data[:,2]), color='k', linestyle='--', linewidth='3', zorder=3, label='GPAW Casida')
        # Manual
        ax.plot([e,e], [0,f[i]], color='C0', linestyle='-', linewidth='2', label='Manual Casida')
    else:
        # GPAW 
        d_e = disc_data[i,1]
        I_e = disc_data[i,2]
        ax.plot([d_e, d_e], [0, I_e] / max(disc_data[:,2]), color='k', linestyle='--', zorder=3, linewidth='2')
        # Manual
        ax.plot([e,e], [0,f[i]], color='C0', linestyle='-', linewidth='2')
ax.set_xlabel('Energy, (eV)')
ax.set_ylabel(r'Photoabsorption spectrum, arb. units')
ax.grid()
plt.tight_layout()
ax.legend(loc='best')
plt.savefig('task2_discrete_spectra.png')

# Finally, convolute my spectrum with a Gaussian and compare with task 1
task1 = np.loadtxt('../task1/spectrum_w0.06.dat', skiprows=4)

energy = np.linspace(0, 7, 500)
f_fold = h.fold(x_t=energy, x_i=eps, y_i=f, width=0.06)
fig, ax = plt.subplots(figsize=(8,6))
ax.plot(task1[:,0], task1[:,1] / max(task1[:,1]), color='k', linestyle='--', linewidth='2', zorder=3, label='GPAW Casida (T.1)') # Normalize both spectra since arbitrary units
ax.plot(energy, f_fold / max(f_fold), color='C0', linestyle='-', linewidth='2', label='Manual Casida') 
ax.set_xlabel('Energy, (eV)')
ax.set_ylabel(r'Photoabsorption spectrum, arb. units')
ax.grid()
plt.tight_layout()
ax.legend(loc='best')
plt.savefig('task2_folded_spectra.png')
plt.show()

