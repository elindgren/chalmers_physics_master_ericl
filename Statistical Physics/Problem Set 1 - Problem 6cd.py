# Imports
import numpy as np
import matplotlib.pyplot as plt

# Set the seed
#np.random.seed(1)

# Number of ensembles
n_ensemble = 100

# Number of steps
n = 1000

# Number of dimensions
ndim = 2

# Pre-allocate vector to save radii in
radii = np.zeros((n_ensemble))

# Iterate over ensembles and calculate ensemble average of radii
for ensemble in range(n_ensemble):
    pass



# Pre-allocate position-vector
R = np.zeros((n,2))
distances = np.zeros((n,1))

# Go through all steps and take a step in some random direction
n_final = 0  # Index at which the chain terminates
last_coord = 0  # From the beginnin, the chain can step in any direction
last_sign = 0
for i in range(2, n+1):
    # Note that the loop starts at 2 - this means that our walk always starts at the origin
    coord = np.random.randint(0,ndim)  # Uniformly choose a coordinate to change
    sign = np.random.randint(low=1,high=2+1)  # Draw a random value between 1 and 2 - will determine if + or - 
    if coord == last_coord and i>2:
        # If its not the first iteration, if the coordinate changed happens to be the same iteration,
        # set the sign to positive so it doesn't go back on itself
        sign = last_sign
    delta = np.zeros((ndim))
    delta[coord] += 1*(-1)**sign
    proposed_position = R[i-1,:] + delta
    # Check if we've already been to this position - if so, terminate
    if np.any((R[:] == proposed_position).all(1)):
        print("Proposed position already visited! Stopping...")
        n_final=i-1  # Last iteration was the final step
        break
    else:
        R[i-1,:] = proposed_position
        last_coord = coord  # The coordinate direction in which a step was taken - can't go opposite to this
        last_sign = sign

    # Calculate the step length for this iteration - from previous step TO THIS iteration
    distances[i-1] = np.linalg.norm(R[i-1,:] - R[i-2,:], ord=2)
    # Set this as our next position
    if not i == n:
        # Don't update next position if currently last point in chain
        R[i,:] = R[i-1,:]
    else:
        n_final = i-1

# Calculate the length of the chain
length=distances.sum()
print(f'Length of chain is {length} units')  # Should be n-1 in units of the step size. 
# Calculate distance from origin:
radius = np.linalg.norm(R[0,:] - R[n_final,:], ord=2)
print(f'Radius of chain is {radius:.2f} units')


# Plot the chain
fig, ax = plt.subplots(figsize=(10,6))
ax.plot(R[:n_final+1,0], R[:n_final+1,1], label=f'Chain length={length:.2f} units')
ax.scatter(R[n_final,0], R[n_final,1], marker='*', c='r', s=100, label=f'Final point')
ax.plot([R[0,0], R[n_final,0]], [R[0,1], R[n_final,1]], label=f'Radius={radius:.2f} units')
ax.legend(loc='best')
plt.show()

# Write a code that simulates a self-avoiding random walk in 2D, and higher dimensions.

