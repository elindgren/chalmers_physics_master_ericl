import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import seaborn as sns



# Set plot params
plt.rc('font', size=14)          # controls default text sizes
plt.rc('axes', titlesize=14)     # fontsize of the axes title
plt.rc('axes', labelsize=14)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=14)    # fontsize of the tick labels
plt.rc('ytick', labelsize=14)    # fontsize of the tick labels
plt.rc('legend', fontsize=14)    # legend fontsize


def free_energy_line(nA, b):
    ''' Calculate the value of the free energy for a given nA and b(eta).'''
    return -0.5*( nA**2 + 0.5*(1-nA)**2 ) + b*( nA*np.log(nA) + (1-nA)*np.log(0.5 * (1-nA)) )


def free_energy_surface(nBs, nCs, b):
    function_map = np.zeros((len(nBs), len(nCs)))
    for i, nB in enumerate(nBs):
        for j, nC in enumerate(nCs):
            if nB+nC <= 1:
                function_map[i,j] = -0.5*( (1-nB-nC)**2 + nB**2 + nC**2 ) + b*( (1-nB-nC)*np.log(1-nB-nC) + nB*np.log(nB) + nC*np.log(nC) ) 
            else:
                function_map[i,j] = np.nan
    return function_map

# Define plot range
nA = np.linspace(0.00001,0.99999,100)
gammas = np.linspace(0,1,10000)

nB = np.linspace(0.00001,0.99999,100)
nC = np.linspace(0.00001,0.99999,100)
Rem = 1-3/2*(nB+nC)
Imm = np.sqrt(3)/2*(nB+nC)

REM, IMM = np.meshgrid(Rem, Imm)

NB, NC = np.meshgrid(nB, nC)
gamma = 1
F = free_energy_surface(nB, nC, gamma)

# Contour plot
fig = plt.figure(figsize=(10,6))
ax = fig.gca(projection='3d')
surf = ax.plot_surface(REM, IMM, F, linewidth=0, antialiased=True)
ax.set_title(rf'Free energy surface, $\gamma=${gamma:.2f}')
ax.set_xlabel(r'Re$(m)$')
ax.set_ylabel(r'Im$(m)$')
ax.set_zlabel(r'Reduced free energy, $\hat{F}$')
plt.tight_layout()
plt.savefig(f'surface_gamma={gamma}.png')
# Find minimum


# Find Tc
fig, ax = plt.subplots(figsize=(10,6))
min_nAs = []

for gamma in gammas:
    # Find minimum for each such and plot nA
    F = free_energy_line(nA, gamma)
    min_idx = np.where(F == np.amin(F))[0]
    min_nAs.append(nA[min_idx])

# Find critical value of gamma
crit_idx = [idx for idx, i in enumerate(min_nAs) if i < 0.4][0]


m = (3*np.array(min_nAs)-1)/2
ax.plot(gammas, m, linewidth=3, color='C2')
ax.scatter(gammas[crit_idx], m[crit_idx], c='C4', s=100, marker='*')
ax.scatter(gammas[crit_idx-1], m[crit_idx-1], c='C4', s=100, marker='*')
ax.annotate(rf'$\Delta m \approx${m[crit_idx-1][0]:.2f}', xy=(gammas[crit_idx-1], m[crit_idx-1]),  xycoords='data',
            xytext=(0.6, 0.6), textcoords='axes fraction',
            arrowprops=dict(facecolor='black', shrink=0.05, width=2),
            horizontalalignment='right', verticalalignment='top',
            )
ax.annotate(rf'$\gamma_c \approx${gammas[crit_idx]:.3f}', xy=(gammas[crit_idx], m[crit_idx]),  xycoords='data',
            xytext=(0.6, 0.2), textcoords='axes fraction',
            arrowprops=dict(facecolor='black', shrink=0.05, width=2),
            horizontalalignment='right', verticalalignment='top',
            )
ax.grid()
ax.set_xlabel(r'$\gamma = \frac{k_B T}{Jz}$')
ax.set_ylabel(r'Magnetization $m_0$')
plt.tight_layout()
plt.savefig(f'transition.png')




# Define plots

# ax.contour3D(N, B, F, 50, cmap='binary')
# ax.plot3D(nA, b, free_energy(nA, b))

plt.show()


