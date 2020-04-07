# System packages
import os
import sys


'''
    A python script for calculating the Young's modulus of the MgO example
    from class. It sets up the necassary LAMMPS files, launches the calculations
    and extracts the information from the results,'.

    Developed by Eric Lindgren, SNUID ericlin 
    April 2020
'''


# Material parameters
a = 4.21079  # Å, relaxed lattice parameter
d = 0.001  # Stretch lattice parameters 0.1% at a time

# Setup
user = 'c2020spring'
axis = 0 # Axis 0 for x, 1 for y, 2 for z
N = 2

p = []
l = []

# Iterate over the various strains
for i in range(-N,N+1):
    ad = a*(1+d*i)  # Strained lattice parameter
    filename = f'MgO-{(1+d*i):.3f}'

    # Setup MgO structure file
    with open('res/MgO.txt', 'r') as file:
        MgO_text = file.readlines()
    row1 = MgO_text[0].split(' ') 
    row1[axis] = f'{ad}' # Modify the lattice parameter for current axis
    MgO_text[0] = ' '.join(row1)
    with open(f'out/{filename}.txt', 'w') as file:
        for row in MgO_text:
            file.write(row)

    # Create data file
    datafile = f'data.{filename}-sc1'
    # os.system(f'/home/{user}/share/lammps-data.exe out/{filename}.txt 1 1 1 > out/f{datafile}')

    # Edit input file 
    with open('res/input.MgO-opt', 'r') as file:
        MgO_text = file.readlines()
    row1 = MgO_text[0].split(' ') 
    row1[-1] = f'{datafile}' # Modify the lattice parameter for current axis
    MgO_text[0] = ' '.join(row1)
    with open(f'out/input.MgO-opt-modified', 'w') as file:
        for row in MgO_text:
            file.write(row)

    # Launch calculation with mpirun
    # os.system(f'mpirun -np 2 /home/bin/lmp_mpich < out/input.MgO-opt-modified > out/output.MgO')

    # Extract final pressure and lattice vectors
    with open(f'out/output.MgO-opt-0.998', 'r') as file:
        output_lines = file.readlines()

    for k, row in enumerate(output_lines):
        if('Loop time' in row):
            info = output_lines[k-1].split('    ')  # Final structure information row in output file
            ps_str = info[2:-3]
            ls_str = info[-3:]
    
    p_ij = []
    l_ij = []
    for j in range(len(ps_str)):
        p_tmp = ps_str[j].split(' ')[0]
        p_ij.append( float(p_tmp) )
        l_ij.append( float(ls_str[j]) )
    
    p.append(p_ij)
    l.append(l_ij)

# Calculate Young's modulus and Poisson ratio
c_ii = (p[0][axis]-p[-1][axis]) / (-2*d - 2*d)
print(c_ii)
print(p)
print(l)