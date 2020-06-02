# Internal imports
from cycler import cycler

# External imports
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import interpolate

# Set plot params
plt.rc('font', size=16)          # controls default text sizes
plt.rc('axes', titlesize=16)     # fontsize of the axes title
plt.rc('axes', labelsize=16)     # fontsize of the x and y labels
plt.rc('xtick', labelsize=16)    # fontsize of the tick labels
plt.rc('ytick', labelsize=16)    # fontsize of the tick labels
plt.rc('legend', fontsize=14)    # legend fontsize
line_cycler = cycler('linestyle',['-','--',':','-.', '-', '--']) + cycler('color', ['b', 'r', 'g', 'k', 'c', 'y'])
plt.rc('axes', prop_cycle=line_cycler)


# Setup
high_cut = 5000 # eV
low_cut = 500 # eV

####* 1. Read digitized csv data
df = pd.read_csv('digitized_pka_data.csv', names=['E', 'CDF'])
E = df['E']
logE = np.log10(df['E'])


####* 2. Spline together to generate N equidistant points
N = 100
tck = interpolate.splrep(logE, df['CDF'], s=0.01)
logenew = np.linspace(0, np.max(logE), N)
enew = np.linspace(1, np.max(E), N)
cdf = interpolate.splev(logenew, tck, der=0)
dx = logenew[1]-logenew[0]

pdf = interpolate.splev(logenew, tck, der=1)



print(f'Integral of PDF in log space: {np.trapz(pdf, logenew):.3f}')
print(f'Integral of PDF in real space: {np.trapz(pdf, enew):.3f}')

####* Monte Carlo Sample CDF
#! CDF is linearly spaced in log10 space - if it was linearly sampled, we would basically never sample the interesting parts of the PDF
Np = 10000
def sample_cdf(logE, cdf, N, seed):
    np.random.seed(seed)
    samples = np.zeros(N)
    rvs = np.random.uniform(low=0, high=1, size=N)
    for i in range(N):
        rv = rvs[i]
        j = 0
        while(rv > cdf[j] and j<len(cdf)-1):
            j += 1
        samples[i] = logE[j]
    return samples
samples = sample_cdf(logenew, cdf, Np, seed=10)

# Write samples with energy less than cutoff to file 
Ns = 50 # Number of samples to write
with open('mc_pka_energies', 'w') as f:
    N_written = 0
    i = 0
    while N_written < Ns:
        log10E = samples[i]
        E_s = 10**log10E
        if E_s < high_cut and E_s > low_cut:
            f.write(f'{E_s:.6f}\n')
            N_written += 1
        i += 1


# Plot
fig, ax = plt.subplots(3,1,figsize=(10,8))
ax[0].plot(logE, df['CDF'], color='k', alpha=0.3, linewidth=2)
ax[0].scatter(logE, df['CDF'], color='k', alpha=0.3, marker='o', s=12, label='Raw data')
ax[0].plot(logenew, cdf, color='b', alpha=0.3, linewidth=2)
ax[0].scatter(logenew, cdf, marker='o', s=12, label='Interpolated data')
ax[0].axvline(np.log10(high_cut), color='r', linestyle='--', label='Cutoff energies')
ax[0].axvline(np.log10(low_cut), color='r', linestyle='--')
ax[0].set_xlabel('log10 of W PKA energy')
ax[0].set_ylabel('CDF')
ax[0].grid()
ax[0].legend(loc='best')


ax[1].plot(logenew, pdf, color='b', alpha=0.3, linewidth=2)
ax[1].scatter(logenew, pdf, marker='o', s=12, label='Interpolated PDF')
ax[1].hist(samples, density=True, bins=100, color='c', label='Samples of PDF')
ax[1].axvline(np.log10(high_cut), color='r', linestyle='--', label='Cutoff energies')
ax[1].axvline(np.log10(low_cut), color='r', linestyle='--')
ax[1].set_xlabel('log10 of W PKA energy')
ax[1].set_ylabel('PDF')
ax[1].grid()
ax[1].legend(loc='best')

normpdf = pdf / np.trapz(pdf, enew)
ax[2].plot(enew, normpdf, color='b', alpha=0.3, linewidth=2)
ax[2].scatter(enew, normpdf, marker='o', s=12)
ax[2].hist(10**samples, density=True, bins=100)
ax[2].set_xlabel('W PKA energy, (eV)')
ax[2].set_ylabel('PDF')
ax[2].grid()

plt.tight_layout()
plt.show()