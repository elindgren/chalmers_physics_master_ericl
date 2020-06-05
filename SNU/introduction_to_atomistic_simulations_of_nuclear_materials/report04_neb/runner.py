# Internal imports
import os
import time

'''
    This script runs the NEB method for finding the migration barrier in 
    4x4x4 supercells of Fe with different 110 dumbbell SIA configurations. 
    The energy calculations then done by analyzer.py.

    Developed by Eric Lindgren, SNU ID ericlin
    May 2020
'''

# Settings
n_img = 21

print(f'Starting NEB calculations for {n_img} images')
s_time = time.time()
for d in os.listdir('./'):
    if os.path.isdir(f'./{d}'):
        d_time = time.time()
        print(f'\tDirectory: {d}', end='')
        os.chdir(f'./{d}')
        # Optimize Fe structures
        os.system('mpirun -np 2 /home/bin/lmp_mpich < input.Fe-relax-ini > tmp')
        os.system('mpirun -np 2 /home/bin/lmp_mpich < input.Fe-relax-fin > tmp')
        # Put fin image in coord.fin
        with open('dump.fin', 'r') as f:
            coords = f.readlines()
        with open('coord.fin', 'w') as f:
            for i, line in enumerate(coords):
                if 'ITEM: NUMBER OF ATOMS' in line:
                    n_atoms = coords[i+1]
            f.write(n_atoms)
            f.writelines(coords[9:])
        # Change the number of images
        with open('../input.neb', 'r') as f:
            inp = f.readlines()
        with open('./input.neb', 'w') as f:
            for line in inp:
                if 'variable num_neb_image equal' in line:
                    f.write(f'variable num_neb_image equal {n_img}\n')
                else:
                    f.write(line)
        # Launch the calculation
        os.system(f'mpirun -np {n_img} /home/bin/lmp_mpich -partition {n_img}x1 -in input.neb > output.neb')
        # Extract reaction coordinates and write to transtion.txt
        with open('output.neb', 'r') as f:
            outp = f.readlines()
        trans = list(filter( None, outp[-1].rstrip().split(' ') ))[9:]
        if d == 'casea':
            mode = 'w'
        else:
            mode = 'a'
        with open('../transitions.txt', mode) as f:
            f.write(f'{d}\n')
            for i in range(n_img):
                f.write(f'{trans[2*i]} {trans[2*i+1]}\n')
        # Return to parent directory
        print(f' --- Finished in {(time.time()-d_time):.2f} s')
        os.chdir(f'../')

print('********')
print(f'Total calculation time: {(time.time()-s_time):.2f} s')