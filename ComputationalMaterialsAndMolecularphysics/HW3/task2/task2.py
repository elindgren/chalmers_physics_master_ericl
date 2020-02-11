# Built-in packages
import os.path

# External packages
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from ase.io.trajectory import Trajectory
from tqdm import tqdm


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


distances = [] # Holds all distances from Na+ for all snapshots 
# Load the trajectory from task 1
traj = Trajectory('../task1/mdTask1.traj')

if not os.path.exists('distances.npy'):
    print('\t----Createing distances vector----')
    for snapshot in tqdm(traj):
        # Get indices for Na+ ion and oxygen atoms
        O_idx = [O.index for O in snapshot if O.symbol=='O']
        Na_idx = [Na.index for Na in snapshot if Na.symbol=='Na']

        # Calculate their distances 
        distances.extend(snapshot.get_distances(Na_idx, indices=O_idx, mic=True))  # Enable minimum image convention to check pbc
    distances = np.array(distances)
    print('\t----Saving distances.npy----')
    np.save('distances.npy', distances)
else:
    print("\t----Load distances.npy----")
    distances = np.load('distances.npy')



# Plot RDF as histogram
hist, r = np.histogram(distances, bins=100)  # Return the bin edges and use the as r vector for normalization
print(r)

# Normalize the obtained distribution function RDF(r) with the number density times the volume of the spherical shell at radius (r)
# This makes sure that the integrated RDF is 1.
rho = 24+1  # Number density of the system - 24 O + 1 Na
hist /= (rho*4*np.pi*r**2)

# fig, ax = plt.subplots(figsize=(8,6))
# ax.hist(distances, bins=100)
# plt.show()