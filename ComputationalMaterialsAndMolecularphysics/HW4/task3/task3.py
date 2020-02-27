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
d = np.load('../task2/dumpTask1.npz')
K = d['K_pp'] # K-matrix
omega = d['ediff_p'] # KS eigenvalue difference
n = d['fdiff_p']  # Occupation difference
mux_p = d['mux_p']
muy_p = d['muy_p']
muz_p = d['muz_p']
mu_p = [mux_p, muy_p, muz_p]

# Construct the Omega matrix
Omega = np.diag(v=omega**2, k=0) # First add diagonal
# Nothing more, since K=0

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
        f_ia = 2 * np.abs( np.sum( mua_p * np.sqrt(n*omega) * F[:,i] ) )
        f[i] += f_ia
    f[i] /= 3 # Average over all dipole moments

# Compare with the discrete spectrum from GPAW
# Plot
# fig, ax = plt.subplots(figsize=(8,6))
# for i,e in enumerate(eps):
#     if i==0:
#         ax.plot([e,e], [0,f[i]], color='C0', linestyle='-', linewidth='2', label='Casida (manual)')
#     else:
#         ax.plot([e,e], [0,f[i]], color='C0', linestyle='-', linewidth='2')
# ax.set_xlabel('Energy, (eV)')
# ax.set_ylabel(r'Oscillator strength, arb. units')
# ax.grid()
# plt.tight_layout()
# plt.savefig('task3_discrete_kohn_sham_spectra.png')


# Finally, convolute my spectrum with a Gaussian and compare with KS eigenvalue differences
print('**** Kohn-Sham eigenvalue differnces ****')
print(omega)
print('*****************************************')
# TODO: What do these KS transitions mean?
energy = np.linspace(0, 7, 1000)
f_fold = h.fold(x_t=energy, x_i=eps, y_i=f, width=0.06)
fig, ax = plt.subplots(figsize=(8,6))

ax.plot(energy, f_fold / max(f_fold), color='C0', linestyle='-', linewidth='2', label='Manual Casida') 
ax.set_xlabel('Energy, (eV)')
ax.set_ylabel(r'Oscillator strength, arb. units')
ax.grid()
plt.tight_layout()
plt.savefig('task3_kohn_sham_spectra.png')
plt.show()

