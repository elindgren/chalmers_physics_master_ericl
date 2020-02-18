# Built-in packages
import os.path

# External packages
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import medfilt
from ase.io.trajectory import Trajectory
from tqdm import tqdm

# Set plot params
plt.rc('font', size=18)          # controls default text sizes
plt.rc('axes', titlesize=18)     # fontsize of the axes title
plt.rc('axes', labelsize=18)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
plt.rc('ytick', labelsize=18)    # fontsize of the tick labels
plt.rc('legend', fontsize=18)    # legend fontsize

'''
Calculate the first ion solvation shell for my simulated trajectory from task 1 and
from the given trajectory in NaCluster24.traj.
The RDF is calculated from the position of the Na-ion (i.e. Na is at r=0) since we 
are interested in the solvation shell of the ion. 

The relevant species for the partial-RDF is the Na+ ion and the oxygen molecules,
due to the hydrogen atoms being bound to the oxygen. Thus the positions of the
hydrogen atoms needs to be integrated out.

The RDF is given as a histogram over the distances from the Na+ ion to the oxygen 
atoms for all snapshots. 
'''

def load_distances(file, trajectory):
    distances = [] # Holds all distances from Na+ for all snapshots 
    if not os.path.exists(f'{file}_{len(trajectory)}_steps.npy'):
        print('---- Creating distances vector ----')
        for snapshot in tqdm(trajectory):
            # Get indices for Na+ ion and oxygen atoms
            O_idx = [O.index for O in snapshot if O.symbol=='O']
            Na_idx = [Na.index for Na in snapshot if Na.symbol=='Na']

            # Calculate their distances 
            distances.extend(snapshot.get_distances(Na_idx, indices=O_idx, mic=True))  # Enable minimum image convention to check pbc
        distances = np.array(distances)
        print('---- Saving to ' + f'{file}_{len(trajectory)}_steps.npy ----'.rjust(20))
        np.save(f'{file}_{len(trajectory)}_steps.npy', distances)
    else:
        print('---- Load from ' + f'{file}_{len(trajectory)}_steps.npy ----'.rjust(30))
        distances = np.load(f'{file}_{len(trajectory)}_steps.npy')
    return distances


def generate_partial_RDF(distances, n_snapshots):
    n = 200  # Number of points 
    # Get a plot of the RDF from its histogram - it is accurate for sufficiently small bins
    RDF, b = np.histogram(distances, bins=n)  # Return the bin edges and use the as r vector for normalization
    RDF = RDF.astype(float) / n_snapshots  # Get the average occupation in each box over all snapshots
    r = np.array([(b[i-1]+b[i])/2 for i in range(1, len(b))])  # Position is middle point of each bin
    dr = r[1]-r[0]  # The size of the spherical radial shell

    # Normalize the obtained distribution function RDF(r) with the number density times the volume of the spherical shell at radius (r)
    # This would have been the RDF had the particles been uncorrelated (on average an rho particles per unit shell).
    V = traj[0].get_volume()
    rho = 24/V  # Number density of the relevant species of the system - ( 24 O + 1 Na ) / volume of unit cell
    RDF /= (rho*4*np.pi*r**2*dr)

    # Pad RDF and r with zeros for a nice looking plot
    r = np.insert(arr=r, obj=0, values=np.linspace(0,r[0],10))
    RDF = np.insert(arr=RDF, obj=0, values=[0]*10)

    return r, RDF


def solv_shell(r, RDF):
    ##### Find the first minimum #####
    
    # Use a median filtered version of the signal to get good estimates on the extremum indices
    filt_RDF = medfilt(RDF, 9)
    # Find the first and the halfway idx, and search inbetween
    first_max = np.argmax(filt_RDF)
    halfway_idx = int(len(filt_RDF)/2)
    min_idx = first_max + np.argmin(filt_RDF[first_max:halfway_idx])

    # Integrate up to that point - the padding does nothing, since it is zero
    shell_size = np.trapz(y=r[:min_idx]*RDF[:min_idx], x=r[:min_idx])
    print(f'First solvation shell: {shell_size:.4f} Å')  # Dimension Å since integral of dimless RDF (just a histogram)
    # Return first solvation shell
    return first_max, halfway_idx, min_idx, shell_size

print("#### Task 2 - Calculate partial RDF ####")
# Load the trajectory from task 1
traj = Trajectory('../task1/mdTask1.traj')
traj_given = Trajectory('../task1/Na-aimd/NaCluster24.traj')

# Calculate the distances between the relevant species
eq_idx = 1000 # Index for equlibration
eq_idx_g = 1000
distances = load_distances('distances', traj[eq_idx:])
distances_given = load_distances('distances_given', traj_given[eq_idx_g:])

# Generate RDF
r, RDF = generate_partial_RDF(distances, len(traj[eq_idx:]))
r_g, RDF_g = generate_partial_RDF(distances_given, len(traj_given[eq_idx_g:]))

# Calculate the first solvation shell of Na for both trajectories
first_max, halfway_idx, min_idx, _ = solv_shell(r, RDF)
max_idx = [first_max, halfway_idx]
_, _, _, _ = solv_shell(r_g, RDF_g)

# Plot
fig, ax = plt.subplots(figsize=(10,6))
ax.plot(r, RDF, linewidth=2, linestyle='-', alpha=1, label=r'Generated trajectory, $2 \rm\, ps$.')
ax.plot(r, medfilt(RDF, 9), linewidth=2, linestyle='-', label=r'Median filtered RDF, $2 \rm\, ps$.')
ax.plot(r_g, RDF_g, linewidth=2, linestyle='--', alpha=1, label=r'Given trajectory, $7 \rm\, ps$.')
# ax.scatter(r[max_idx], RDF[max_idx], marker='o', s=48, c='k')
# ax.scatter(r[min_idx], RDF[min_idx], marker='s', s=48, c='k')
ax.axhline(1, linestyle='--', c='k')
ax.axvline(r[min_idx], linestyle=':', c='k', label=r'First min($g_{NO}$), $2 \rm\, ps$.')
ax.legend(loc='best')
ax.set_xlabel(r'$r$, (Å)')
ax.set_ylabel(r'$g_{\rm Na^+O}$')
ax.set_xlim(0,6)
ax.grid()
plt.tight_layout()
plt.savefig('task2.png')
plt.show()

print("#### Task 2 - " +  "Finished ####".rjust(26))
