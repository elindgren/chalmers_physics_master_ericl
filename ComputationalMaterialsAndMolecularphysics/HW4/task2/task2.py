# External imports
import numpy as np
import matplotlib.pyplot as plt

# Load dumped data from task 1

# Construct the Omega matrix
Omega = np.diag(v=omega, k=0) # First add diagonal
for p in range(len(Omega)):
    for q in range(len(Omega)):
        Omega[p,q] += 2*np.sqrt(n[p]*omega[p]) * K[p,q] * np.sqrt(n[q]*omega[q])


# Obtain egivenalues and eigenvectors from the Omega matrix

# Calculate the oscillator strength

# Compute the discrete spectrum

# Compare with the discrete spectrum from GPAW

# Finally, convolute my spectrum with a Gaussian and compare with task 1