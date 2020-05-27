# Internal imports
import os
import time
import math
import random

# # Comment out these
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
# import numpy as np
# # Set plot params
# plt.rc('font', size=16)          # controls default text sizes
# plt.rc('axes', titlesize=14)     # fontsize of the axes title
# plt.rc('axes', labelsize=14)     # fontsize of the x and y labels
# plt.rc('xtick', labelsize=14)    # fontsize of the tick labels
# plt.rc('ytick', labelsize=14)    # fontsize of the tick labels
# plt.rc('legend', fontsize=14)    # legend fontsize

'''
    This script runs the LAMMPS calculations for the collision cascade study of radiation damage
    in 20x15x15 atom supercell of bcc-Fe. 

    Developed by Eric Lindgren, SNU ID ericlin
    May 2020
'''

def first_quadrant_20():
    ''' Returns 20 directions in the first quadrant of a sphere with readius r = 1.'''
    pi = 3.14159265359
    # Theta - inclination angle, Phi - azimuthal angle
    nt = 4 
    np = 5
    directions = []
    for t in range(nt):
        for p in range(np):
            theta = random.random()*pi/2
            phi = random.random()*pi/2
            x = math.sin(theta)*math.cos(phi)
            y = math.sin(theta)*math.sin(phi)
            z = math.cos(theta)
            directions.append([x,y,z])
    return directions


# LAMMPS controls
random.seed(2)
PKA = 1  # AtomID of PKA
directions = first_quadrant_20()

# Plot directions just for clarity
# directions = np.array(directions)
# fig = plt.figure(figsize=(8,6))
# ax = fig.add_subplot(111, projection='3d')
# for d in directions:
#     ax.scatter(d[0], d[1], d[2], color='k')
#     ax.plot(xs=[d[0], 0], ys=[d[1], 0], zs=[d[2], 0], color='k', linestyle='--', alpha=0.3, linewidth=1)
# ax.set_xlabel('x')
# ax.set_ylabel('y')
# ax.set_zlabel('z')
# plt.tight_layout()
# plt.savefig('directions.png')

start_time = time.time()
for i, d in enumerate(directions):
    d_time = time.time()
    print(f'\t{i+1}/20 - Direction: ({d[0]:.2f}, {d[1]:.2f}, {d[2]:.2f})', end='')
    # Make new folder
    fdr = f'direction{i+1}'
    if not os.path.exists(fdr):
        os.system(f'mkdir {fdr}')
    os.chdir(f'./{fdr}')
    # Modify input script to change PKA atom and direction
    with open('../input.cascade-Fe', 'r') as f:
        inp = f.readlines()
    with open('./input.cascade-Fe-mod', 'w') as f:
        for row in inp:
            if 'variable pkaid1 equal' in row:
                f.write(f'variable pkaid1 equal {PKA}\n')
            elif 'variable pkadx equal' in row:
                f.write(f'variable pkadx equal {d[0]}\n')
            elif 'variable pkady equal' in row:
                f.write(f'variable pkadx equal {d[1]}\n')
            elif 'variable pkadz equal' in row:
                f.write(f'variable pkadz equal {d[2]}\n')
            else:
                f.write(row)
    # Start calculation 
    os.system('mpirun -np 2 /home/bin/lmp_mpich < input.cascade-Fe-mod > output.cascade-Fe-mod')
    # Move back up directory structure
    print(f' --- Finished in {(time.time()-d_time):.2f} s')
    os.chdir('../')
print('********')
print(f'Total calculation time: {(time.time()-start_time):.2f} s')