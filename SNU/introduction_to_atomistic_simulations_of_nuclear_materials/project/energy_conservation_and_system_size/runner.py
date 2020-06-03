# Internal imports
import time
import shutil
import os
import sys

'''
    This script runs energy and system size tests for PKA simulations in W
    with a PKA energy of 5 keV with different timesteps.
    Abort calculation after T picoseconds
'''

dts = [
    0.002,
    0.001,
    0.0005,
    0.0004,
    0.0003,
    0.0002,
    0.0001,
]
xs = [
    0.01,
    0.005,
    0.0025,
    0.0020,
    0.0015,
    0.0010,
    0.0005
]

Ncpu = 8


def main(arg):
    try:
        pkaE = int(arg[0])
    except:
        print(f'Incorrect PKA Energy {arg}')
        sys.exit(1)
    Edir = f'./{pkaE}eV'
    if not os.path.isdir(Edir):
        os.mkdir(Edir) # Create folder
    os.chdir(Edir) # Move to folder

    edirect = f'./energies_{pkaE}eV'
    if os.path.isdir(edirect):
        shutil.rmtree(edirect)
    os.mkdir(edirect) # Make a new energies folder

    # Start script
    print(f'PKA Energy: {pkaE} eV')
    st = time.time()
    for i, dt in enumerate(dts):
        d_time = time.time()
        print(f'\tTime step: {dt}', end='')

        x = xs[i]
        direct = f'./dt={dt}'
        if not os.path.isdir(direct):
            os.mkdir(direct) # Create folder
        os.chdir(direct) # Move to folder
        # Create edited input script
        with open('../../res/input.cascade-energy', 'r') as f:
            inp = f.readlines()
        with open('input.cascade-energy', 'w') as f:
            for line in inp:
                if 'variable pkaene equal' in line:
                    f.write(f'variable pkaene equal {pkaE}\n')
                elif 'variable tmax1 equal' in line:
                    f.write(f'variable tmax1 equal {dt}\n')
                elif 'variable xmax1 equal' in line:
                    f.write(f'variable xmax1 equal {x}\n')
                else:
                    f.write(line)
        resdirect = f'./restart'
        if os.path.isdir(resdirect):
            shutil.rmtree(resdirect)
        os.mkdir(resdirect) # Make a new restart folder - remove old restart files
        os.system(f'mpirun -np {Ncpu} /home/bin/lmp_mpich < input.cascade-energy > output.cascade-energy') # Launch calculation
        print(f' --- Finished in {(time.time()-d_time):.2f} s on {Ncpu} cores')
        os.chdir(f'../') # Return to parent directory
    print('********')
    print(f'Total calculation time: {(time.time()-st):.2f} s on {Ncpu} cores')

if __name__ == "__main__":
   main(sys.argv[1:])