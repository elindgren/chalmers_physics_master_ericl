# Built-in packages
import time

# ASE
from ase.parallel import world

# GPAW
from gpaw import GPAW
from gpaw.lrtddft import LrTDDFT
from gpaw.lrtddft import photoabsorption_spectrum


if world.rank==0:
    print(f'------------   Extracting photoabsorption spectrum using TDDFT   ------------')
start = time.time()

# Load calculator after relaxing empty structure
calc = GPAW('emptyCalc.gpw')
calc.set(
    txt = 'spec.gpaw-out'
)

# Calculate and diagonalize Omega matrix
dE = 6  # Up to 6 eV transitions considered
lr = LrTDDFT(
    calc, 
    xc='LDA', 
    energy_range=dE,
)  # Construct the omega matrix, parallelised over all available cores

lr.write(f'lr_dE={dE}eV.dat.gz')  # Save the tdDFT calculater just in case

lr.diagonalize()
wd = 0.06
photoabsorption_spectrum(
    lr,
    f'spectrum_w{wd}.dat',
    width = wd
)


 # Save results for task 2
lr.write ('TDDFT_Task1.dat')


end = time.time()
if world.rank==0:
    print('-------- Photoabsorption spectrum extracted in: ' + f'{(end-start):.2f} s --------'.rjust(34))
    print('----------------------------------------------------------------------')