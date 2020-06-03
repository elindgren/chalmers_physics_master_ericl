# Internal imports
import time
import os

'''
    This script runs energy and system size tests for PKA simulations in W
    with a PKA energy of 5 keV with different timesteps.
    Abort calculation after T picoseconds
'''

dts = [
    0.002,
    0.001,
    0.0005,
    0.0001
]
xs = [
    0.01,
    0.005,
    0.0025,
    0.0005
]

Ncpu = 8


st = time.time()
for i, dt in enumerate(dts):
    d_time = time.time()
    print(f'\tTime step: {dt}', end='')

    x = xs[i]
    direct = f'./dt={dt}'
    if not os.path.isdir(direct):
        os.mkdir(direct) # Create folder
    os.chdir(f'./dt={dt}') # Move to folder
    # Create edited input script
    with open('../res/input.cascade-energy', 'r') as f:
        inp = f.readlines()
    with open('input.cascade-energy', 'w') as f:
        for line in inp:
            if 'variable tmax1 equal' in line:
                f.write(f'variable tmax1 equal {dt}\n')
            elif 'variable xmax1 equal' in line:
                f.write(f'variable xmax1 equal {x}\n')
            else:
                f.write(line)
    resdirect = f'./restart'
    if not os.path.isdir(resdirect):
        os.mkdir(resdirect) # Make restart folder
    os.system(f'mpirun -np {Ncpu} /home/bin/lmp_mpich < input.cascade-energy > output.cascade-energy') # Launch calculation
    print(f' --- Finished in {(time.time()-d_time):.2f} s on {Ncpu} cores')
    os.chdir(f'../') # Return to parent directory
print('********')
print(f'Total calculation time: {(time.time()-st):.2f} s on {Ncpu} cores')