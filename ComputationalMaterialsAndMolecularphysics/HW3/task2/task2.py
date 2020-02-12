# Built-in packages
import os.path

# External packages
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from ase.io.trajectory import Trajectory
from tqdm import tqdm

# Set plot params
plt.rc('font', size=14)          # controls default text sizes
plt.rc('axes', titlesize=14)     # fontsize of the axes title
plt.rc('axes', labelsize=16)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=14)    # fontsize of the tick labels
plt.rc('ytick', labelsize=14)    # fontsize of the tick labels
plt.rc('legend', fontsize=14)    # legend fontsize

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
    if not os.path.exists(f'{file}.npy'):
        print('---- Creating distances vector ----')
        for snapshot in tqdm(trajectory):
            # Get indices for Na+ ion and oxygen atoms
            O_idx = [O.index for O in snapshot if O.symbol=='O']
            Na_idx = [Na.index for Na in snapshot if Na.symbol=='Na']

            # Calculate their distances 
            distances.extend(snapshot.get_distances(Na_idx, indices=O_idx, mic=True))  # Enable minimum image convention to check pbc
        distances = np.array(distances)
        print('---- Saving to ' + f'{file}npy ----'.rjust(20))
        np.save(f'{file}.npy', distances)
    else:
        print('---- Load from ' + f'{file}.npy ----'.rjust(30))
        distances = np.load(f'{file}.npy')
    return distances


def generate_partial_RDF(distances):
    n = 200  # Number of points 
    # Get a plot of the RDF from its histogram - it is accurate for sufficiently small bins
    RDF, b = np.histogram(distances, bins=n)  # Return the bin edges and use the as r vector for normalization
    RDF = RDF.astype(float)
    r = np.array([(b[i-1]+b[i])/2 for i in range(1, len(b))])  # Position is middle point of each bin
    dr = r[1]-r[0]  # The size of the spherical radial shell

    # Normalize the obtained distribution function RDF(r) with the number density times the volume of the spherical shell at radius (r)
    # This would have been the RDF had the particles been uncorrelated (on average an rho particles per unit shell).
    V = traj[0].get_volume()
    rho = 25/V  # Number density of the relevant species of the system - ( 24 O + 1 Na ) / volume of unit cell
    RDF /= (rho*4*np.pi*r**2*dr)

    # Pad RDF and r with zeros for a nice looking plot
    r = np.insert(arr=r, obj=0, values=np.linspace(0,r[0],10))
    RDF = np.insert(arr=RDF, obj=0, values=[0]*10)

    return r, RDF


print("#### Task 2 - Calculate partial RDF ####")
# Load the trajectory from task 1
traj = Trajectory('../task1/mdTask1.traj')
traj_given = Trajectory('../task1/Na-aimd/NaCluster24.traj')

# Calculate the distances between the relevant species
distances = load_distances('distances', traj)
distances_given = load_distances('distances_given', traj_given)

# Generate RDF
r, RDF = generate_partial_RDF(distances)
r_g, RDF_g = generate_partial_RDF(distances_given)


# Plot
fig, ax = plt.subplots(figsize=(8,6))
ax.plot(r, RDF, linewidth=2, linestyle='-', label=r'Generated trajectory, $2 \rm\, ps$.')
ax.plot(r_g, RDF_g, linewidth=2, linestyle='--', label=r'Given trajectory, $7 \rm\, ps$.')
ax.axhline(1, linestyle=':', c='k')
ax.legend(loc='best')
ax.set_xlabel(r'$r$, (Å)')
ax.set_ylabel(r'$g_{NO}$')
ax.grid()
plt.tight_layout()
plt.show()

print("#### Task 2 - " +  "Finished ####".rjust(26))
