# Imports
import numpy as np
import matplotlib.pyplot as plt

# Set the seed
np.random.seed(1)

# Number of steps
n = 1000

# Number of dimensions
ndim = 2

# Pre-allocate position-vector
R = np.zeros((n,2))
distances = np.zeros((n,1))

# Go through all steps and take a step in some random direction
for i in range(2, n+1):
    # Note that the loop starts at 2 - this means that our walk always starts at the origin
    coord = np.random.randint(0,ndim)  # Uniformly choose a coordinate to change
    sign = np.random.randint(low=1,high=2+1)  # Draw a random value between 1 and 2 - will determine if + or - 
    R[i-1, coord] += 1*(-1)**sign
    # Calculate the step length for this iteration - from previous step TO THIS iteration
    distances[i-1] = np.linalg.norm(R[i-1,:] - R[i-2,:], ord=2)
    # Set this as our next position
    if not i == n:
        # Don't update next position if currently last point in chain
        R[i,:] = R[i-1,:]
    

# Calculate the length of the chain
length=distances.sum()

print(f'Length of chain is {length} units')  # Should be n-1 in units of the step size. 

# Plot the chain
fig, ax = plt.subplots(figsize=(10,6))
ax.plot(R[:,0], R[:,1])

plt.show()

# Write a code that simulates a self-avoiding random walk in 2D, and higher dimensions.

