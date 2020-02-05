# Built-in packages
import time

# Third-party packages
import numpy as np

from ase import Atoms
from ase.io import read, write, Trajectory
from ase.visualize import view
from ase.units import fs, kB
from ase.md.npt import NPT

from gpaw import GPAW 


# Load atoms object
a = read('thermalizedConfiguration.xyz')
# view(a)  # DEBUG
 
# Define timestep, total simulation length and number of steps
dt = 0.5*fs
t_tot = 2000*fs  # 2 ps
N_steps = int(t_tot/dt)
print(f'------------    MD simulation with GPAW for {N_steps} steps    ------------')
start = time.time()
calc = GPAW(mode='lcao',
            xc='PBE',
            basis='dzp',
            symmetry={'point_group': False},
            charge=1,
            txt='out_mdtask1.txt'
)

a.set_calculator(calc)

dyn = NPT(
    atoms = a,
    timestep=dt,
    temperature=350*kB,
    externalstress=0,
    ttime=20*fs,
    pfactor=None,
    logfile='mdTask1.log'
) # Using the Nosé-Hoover thermostat

trajectory = Trajectory('mdTask1.traj', 'w', a)
dyn.attach(trajectory.write, interval=1) # Write state of the system to trajectory at every timestep
dyn.run(N_steps)
end = time.time()
print('-------- MD simulation finished in: ' + f'{(end-start):.2f} s --------'.rjust(34))
print('----------------------------------------------------------------------')