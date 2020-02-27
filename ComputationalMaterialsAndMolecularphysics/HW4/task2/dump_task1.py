# External imports
import helper as h

from gpaw import GPAW
from gpaw.lrtddft import LrTDDFT

# Retrieve the calculator from task 1
lr = LrTDDFT('../task1/TDDFT_Task1.dat')

# Dump the results from the calculation to file
h.dump_data(lr, fpath='dumpTask1.npz' )

# Also generate the discrete spectrum from GPAW
h.discrete_spectrum(lr, 'GPAW_discrete.dat')







