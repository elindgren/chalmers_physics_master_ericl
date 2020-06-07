# External imports
import pandas as pd
import numpy as np

# Load data
def read_force(file):
    forces = []
    with open(file, 'r') as f:
        dmp = f.readlines()
        for i, line in enumerate(dmp):
            if i>8:
                tmp = [s for s in line.rstrip().split(' ') if not s=='']
                forces.append(float(tmp[-3]))
    return np.array(forces)


delta = 0.005

forces_plus = read_force('dump.opt_plus')
forces_minus = read_force('dump.opt_minus')

difference = (forces_plus-forces_minus) / (2*delta)
print(difference)

