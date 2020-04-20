# System packages
import os
import sys
import time
import statistics as s
import time as t


'''
    A python script for generating Fe.txt files with 
    perturbed lattice parameters.

    Developed by Eric Lindgren, SNUID ericlin 
    April 2020
'''


# Material parameters
a = 2.8553  # Å, relaxed lattice parameter
d = 0.001  # Stretch lattice parameters 0.1% at a time

# Setup
user = 'c2020spring'
axis = 0 # Axis 0 for x, 1 for y, 2 for z
N = 2

# Resources -------- Change here!
mat = 'Fe'  # material
st_txt = 'Fe.txt'  # txt structure file
in_file = 'input.Fe-0K'
s = 1  # Size of side of supercell

start_time = t.time()

# Iterate over the various strains
for i in range(-N,N+1):
    print()
    print(f'-------- Iteration {i+3} --------')
    ad = a*(1+d*i)  # Strained lattice parameter
    filename = f'{mat}-{s}x{s}x{s}-{(1+d*i):.3f}'

    # Setup structure file
    with open(f'res/{st_txt}', 'r') as file:
        txt = file.readlines()
    row1 = txt[0].split(' ')
    row1[axis] = f'{ad}' # Modify the lattice parameter for current axis
    txt[0] = ' '.join(row1)
    with open(f'out/{filename}.txt', 'w') as file:
        for row in txt:
            file.write(row)

    # Create data file
    datafile = f'data.{filename}'
    print('\tSetting up datafile')
    os.system(f'/home/{user}/share/lammps-data.exe out/{filename}.txt {s} {s} {s} > out/{datafile}')

    # Edit input file 
    with open(f'res/{in_file}', 'r') as file:
        input_text = file.readlines()
    for i, row in enumerate(input_text):
        if 'variable dfile' in row:
            drow = row.split(' ')
            r = i
            break
    drow[-1] = f'out/{datafile}\n' # Modify the lattice parameter for current axis
    input_text[r] = ' '.join(drow)
    with open(f'out/{in_file}-modified', 'w') as file:
        for row in input_text:
            file.write(row)
    file.close()
        
    # Launch calculation with mpirun
    print('\tLaunching LAMMPS calculation')
    calc_time = t.time()
    os.system(f'mpirun -np 2 /home/bin/lmp_mpich < out/{in_file}-modified > out/output.{filename}')
    print(f'-------- LAMMPS calculation finished in: {(t.time() - calc_time):2f} s --------')



print(f'#########    All calculations finished in: {(t.time() - start_time):2f} s    #########')