# Internal imports
import time
import shutil
import os
import sys

'''
    This script runs mechanical properties calculations for all the resulting structures in the
    run_cascade_simulations folder.
'''

home_dir = os.getcwd()
# print(home_dir)
for obj in os.listdir('../run_cascade_calculations/'):
    Efold = f'../run_cascade_calculations/{obj}'
    if os.path.isdir(Efold) and 'eV' in Efold and not 'not_complete' in Efold:
        E = Efold.split('/')[-1].split('eV')[0]
        print(f'Energy: {E} eV')
        # Remove old mehcanical results file
        resFile = f'{E}_eV_mechanical.results'
        if os.path.isfile(resFile):
            os.remove(resFile)
        os.chdir(Efold) 
        for direction in os.listdir('.'):
            if os.path.isdir(direction):
                d_time = time.time()
                print(f'\tDirection: {direction}', end='')
                os.chdir(direction) # Move to directory and launch calculation
                os.system(f'cp {home_dir}/in.elastic-rev .')
                os.system(f'cp {home_dir}/init.mod .')
                os.system(f'cp {home_dir}/displace.mod .')
                os.system(f'cp {home_dir}/potential.mod .')
                # Modify data.after_event to be triclinic
                with open('data.after_event', 'r') as orig:
                    with open('data.after_event-mod', 'w') as mod:
                        lines = orig.readlines()
                        for line in lines:
                            mod.write(line)
                            if 'zlo zhi' in line:
                                mod.write('0.0 0.0 0.0 xy xz yz\n')
                os.system('mpirun -np 8 /home/bin/lmp_mpich < in.elastic-rev > out.elastic')
                # Read results and write them to file
                bulk = 0
                shear1 = 0
                shear2 = 0
                poisson = 0
                with open('out.elastic', 'r') as e:
                    lines = e.readlines()
                    for line in lines:
                        if 'Bulk Modulus' in line:
                            val = line.rstrip().split(' ')[-2]
                            bulk = float(line.rstrip().split(' ')[-2])
                        elif 'Shear Modulus 1' in line:
                            shear1 = float(line.rstrip().split(' ')[-2])
                        elif 'Shear Modulus 2' in line:
                            shear2 = float(line.rstrip().split(' ')[-2])
                        elif 'Poisson Ratio' in line:
                            poisson = float(line.rstrip().split(' ')[-1])
                # Write results to file
                with open(f'{home_dir}/{resFile}', 'a') as m:
                    d = ''.join(direction.split(', '))
                    m.write(f'{d},{bulk},{shear1},{shear2},{poisson}\n')
                print(f' --- Finished in {(time.time()-d_time):.2f} s')
                os.chdir('..') # Move back to energy dir
                

        