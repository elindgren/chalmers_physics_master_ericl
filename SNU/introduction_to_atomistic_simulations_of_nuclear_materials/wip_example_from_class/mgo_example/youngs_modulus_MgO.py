# System packages
import os
import sys
import time
import statistics as s


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
    print()
    print(f'-------- Iteration {i+3} --------')
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
    datafile = f'data.{filename}'
    print('Setting up datafile')
    os.system(f'/home/{user}/share/lammps-data.exe out/{filename}.txt 1 1 1 > out/{datafile}')

    # Edit input file 
    with open('res/input.MgO-opt', 'r') as file:
        MgO_text = file.readlines()
    row1 = MgO_text[0].split(' ') 
    row1[-1] = f'out/{datafile}' # Modify the lattice parameter for current axis
    MgO_text[0] = ' '.join(row1)
    with open(f'out/input.MgO-opt-modified', 'w') as file:
        for row in MgO_text:
            file.write(row)
    file.close()
        
    # Launch calculation with mpirun
    print('Launching LAMMPS calculation')
    os.system(f'mpirun -np 2 /home/bin/lmp_mpich < out/input.MgO-opt-modified > out/output.{filename}')

    # Extract final pressure and lattice vectors
    with open(f'out/output.{filename}', 'r') as file:
        output_lines = file.readlines()

    for k, row in enumerate(output_lines):
        if('Loop time' in row):
            info = output_lines[k-1].split(' ')  # Final structure information row in output file
            info = list(filter(None, info))
            ps_str = info[1:4]
            ls_str = info[-4:-1]
    
    p_ij = []
    l_ij = []
    for j in range(len(ps_str)):
        p_tmp = ps_str[j].split(' ')[0]
        p_ij.append( float(p_tmp) )
        l_ij.append( float(ls_str[j]) )
    
    p.append(p_ij)
    l.append(l_ij)

# Calculate Young's modulus and Poisson ratio
print('***** Calculating parameters of interest *****')
e_11 = []
e_22 = []
for i in range(len(l)):
    e_11.append( l[i][axis]/l[2][axis] - 1 )  # Strain direction
    e_22.append( l[i][1]/l[2][1] - 1 )  # Other direction
v = -(e_22[0]-e_22[-1]) / (e_11[0] - e_11[-1])
print(f'Poisson ratio: {v:4f}')

c_ii = (p[0][axis]-p[-1][axis]) / (e_11[0] - e_11[-1])
print(f'Youngs Modulus: {(c_ii/1e4):4f} GPa')

