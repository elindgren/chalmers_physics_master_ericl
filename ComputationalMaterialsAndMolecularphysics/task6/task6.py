# Internal imports
import os.path as p

# External imports
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

# ASE
from ase import Atoms
from ase.db import connect
from ase.calculators.eam import EAM
from ase.build import bulk
from ase.vibrations import Vibrations
from ase.visualize import view


''' 
Perform GPAW DFT calculation for the electronic structure of 
bulk Al. First find converge number of k-points. Then, converge
density self-consistently.
'''

# DBs
bulkDB = connect('./bulk.db', append=False)  # DB for vibration spectrum

# Using optimal lattice parameter from task 2
a = 4.043  # A
atoms = bulk('Al', 'fcc', a)


