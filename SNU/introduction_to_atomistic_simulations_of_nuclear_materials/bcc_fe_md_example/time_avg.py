# Internal imports

# External imports
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

'''
    Calculates and plots the time averages of various quantities for the 
    BCC Fe example witha 6x6x6 supercell run on lecture 7 (4/9).
'''

N = 6*6*6  # nbr of atoms
col = 6  # Which column to take time average of 

data = []
steps = []
with open('output-modified.Fe-MD-0300K') as file:
    header = file.readline().split(' ')
    label = header[col]
    for row in file:
        split = list(filter(None, row.split(' ')))
        steps.append(float(split[0]))
        data.append(float(split[col]))

data = np.array(data)
steps = np.array(steps)

# Calculate average per atom
avg = np.mean( data )
print(f'Average {label}: {avg:.4f} - Average {label} per atom: {(avg/N):.4f}')
# plot
fig, ax = plt.subplots(figsize=(8,6))
ax.plot(steps, data)
ax.set_ylabel(label)
plt.show()