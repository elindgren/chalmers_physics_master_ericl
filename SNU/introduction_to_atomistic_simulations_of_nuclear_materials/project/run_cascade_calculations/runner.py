# Internal imports
import time
import shutil
import os
import sys

'''
    This script runs collision cascade computations for different combinations of energy and angle.
    Tunable parameters:
        Ncpu: Number of CPU cores to use
        dt: maximum dynamic timestep as set for dt/reset.
        x: maximum distance for dt/reset.
        tmax: Total simulation time, after which cascade simulation is aborted. Note that first 0.2 ps is warmup. 
        pkaID: AtomID for chosen PKA atom. 
'''

Ncpu = 8
dt = 0.001
x = 0.1
pkaID = 152449
tmax = 8 # ps, total simulation time. Note that first 0.2 ps is warmup. 


# Load sampled energies
pkaEnergies = []
with open('mc_pka_energies', 'r') as f:
    for line in f:
        pkaEnergies.append(float(line.rstrip()))

# Load sampled angles
pka_angle = [
    [1,0,0],
    [1,1,0],
    [1,1,1],
    [3,1,0],
    [3,2,0],
    [3,3,1],
    [3,3,2],
    [3,1,1],
    [3,2,2],
    [3,2,1]
]

st = time.time()

with open('runner.out', 'w') as ro:
    ro.write('----- W COLLISION CASCADE SIMULATION -----\n')
    ro.write('\n')

for i, pkaE in enumerate(pkaEnergies):
    # Start script
    e_time = time.time()
    with open('runner.out', 'a') as ro:
        ro.write(f'---- ({i+1}/{len(pkaEnergies)}) PKA Energy: {pkaE:.3f} eV '.ljust(91, '-') + '\n')
    Edir = f'./{pkaE:.3f}eV'
    if not os.path.isdir(Edir):
        os.mkdir(Edir) # Create folder
    os.chdir(Edir) # Move to folder

    for j, ang in enumerate(pka_angle):
        a_time = time.time()
        ang_str = f'[{ang[0]:.2f},{ang[1]:.2f},{ang[2]:.2f}]'
        with open('../runner.out', 'a') as ro:
            ro.write(f'\t ({j+1}/{len(pka_angle)}) '.ljust(12, '-')  + '--- PKA angle: ' + f'{ang_str}'.rjust(25))
        angDir = f'./{ang}'
        if not os.path.isdir(angDir):
            os.mkdir(angDir) # Create folder
        os.chdir(angDir) # Move to folder
        
        # Create edited input script
        with open('../../res/input.cascade-energy', 'r') as f:
            inp = f.readlines()
        with open(f'input.cascade-energy', 'w') as f:
            for line in inp:
                if 'variable pkaid1 equal' in line:
                        f.write(f'variable pkaid1 equal {pkaID}\n')
                elif 'variable pkaene equal' in line:
                    f.write(f'variable pkaene equal {pkaE:.3f}\n')
                elif 'variable pkadx equal' in line:
                    f.write(f'variable pkadx equal {ang[0]:.3f}\n')
                elif 'variable pkady equal' in line:
                    f.write(f'variable pkady equal {ang[1]:.3f}\n')
                elif 'variable pkadz equal' in line:
                    f.write(f'variable pkadz equal {ang[2]:.3f}\n')
                elif 'variable tmax1 equal' in line:
                        f.write(f'variable tmax1 equal {dt}\n')
                elif 'variable xmax1 equal' in line:
                    f.write(f'variable xmax1 equal {x}\n')
                elif 'variable tmax equal' in line:
                    f.write(f'variable tmax equal {tmax}\n')
                else:
                    f.write(line)
        resdirect = f'./restart'
        if os.path.isdir(resdirect):
            shutil.rmtree(resdirect)
        os.mkdir(resdirect) # Make a new restart folder - remove old restart files
        os.system(f'mpirun -np {Ncpu} /home/bin/lmp_mpich < input.cascade-energy > output.cascade-energy') # Launch calculation
        with open('../../runner.out', 'a') as ro:
            ro.write(f' --- Finished in' +  f'{(time.time()-a_time):.2f} s -'.rjust(20) + '\n')
        os.chdir(f'../') # Return to parent directory
    with open('../runner.out', 'a') as ro:
        ro.write(f'---- Finished in {(time.time()-e_time):.2f} s ----'.ljust(91,'-') + '\n')
        ro.write('\n')
    os.chdir(f'../') # Return to parent directory


with open('runner.out', 'a') as ro:
    ro.write('\n')
    ro.write(''.rjust(92,'*') + '\n')
    ro.write(f'****           Total calculation time: ' + f'{(time.time()-st):.2f} s on {Ncpu} cores        *****'.rjust(53) + '\n')
    ro.write(f'****           PKA ID: {pkaID}           Timestep: ' + f'{dt} ps           xstep: {x} A '.ljust(37) + '*****' + '\n')
    ro.write(''.rjust(92,'*') + '\n')