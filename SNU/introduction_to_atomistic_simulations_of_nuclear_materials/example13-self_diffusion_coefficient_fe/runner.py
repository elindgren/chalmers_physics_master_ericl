# Internal imports
import os
import time

'''
    Calculates MSD for BCC Fe at 600 K, both for a perfect crystal and for
    a SIA 111 defect.
'''


# Define cases
cases = ['perfect', 'SIA']
s_time = time.time()
for case in cases:
    c_time = time.time()
    print(f'Performing calculation for: {case}')
    ####* Relax structure

    # Modify relax input script
    restag = f'Fe-relax-{case}'
    dataf = f'data.Fe-sc4-{case}'
    with open('res/input.Fe-relax', 'r') as re:
        inp = re.readlines()
        with open(f'input.{restag}', 'w') as reinp:
            for line in inp:
                if 'variable dfile string' in line:
                    reinp.write(f'variable dfile string res/{dataf}\n')
                else: 
                    reinp.write(line)
    os.system(f'mpirun -np 2 /home/bin/lmp_mpich < input.{restag} > output.{restag}')
    ####* Run MSD 

    # Modify MSD input script
    msdtag = f'Fe-600K-npt-msd-{case}'
    with open('res/input.Fe-600K-npt-msd', 'r') as msd:
        inp = msd.readlines()
        with open(f'input.{msdtag}', 'w') as msdinp:
            for line in inp:
                if 'variable dfile string' in line:
                    msdinp.write(f'variable dfile string res/{dataf}-gopted\n')
                elif 'variable projectname string' in line:
                    msdinp.write(f'variable projectname string sc4-{case}\n')

                else: 
                    msdinp.write(line)
    os.system(f'mpirun -np 2 /home/bin/lmp_mpich < input.{msdtag} > output.{msdtag}')
    print(f'\t --- Finished in {(time.time()-c_time):.2f} s')

print('********')
print(f'Total calculation time: {(time.time()-s_time):.2f} s ')