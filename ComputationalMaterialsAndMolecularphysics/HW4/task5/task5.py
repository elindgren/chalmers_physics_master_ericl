# Built-in packages
import time

# GPAW
from gpaw import GPAW
from gpaw.lrtddft import LrTDDFT
from gpaw.lrtddft import photoabsorption_spectrum


print(f'------------   Extracting shortened photoabsorption spectrum   ------------')

start = time.time()

# Import LrTDDFT results from Task 1    
lr = LrTDDFT('../task1/TDDFT_Task1.dat')

lr.diagonalize(energy_range=4) # Only include up to 4 eV

# Generate spectrum and save it
wd = 0.06
photoabsorption_spectrum(
    lr,
    f'spectrum_w{wd}.dat',
    width = wd
)

# Extract all information about all transitions
print('** LrTDDFT.analyze() output **')
lr.analyze()
print('*******************************')


end = time.time()
print('-------- Photoabsorption spectrum extracted in: ' + f'{(end-start):.2f} s --------'.rjust(34))
print('----------------------------------------------------------------------')