# Internal imports
import os 
import pickle



sc = [3,4,5] # Side of supercell
Ed = {}


for s in sc:
    Ed[f'{s}'] = {}
    d = f'./{s}x{s}x{s}'
    if not os.path.isdir(d):
        os.mkdir(d)
    # Change dir
    os.chdir(d)
    # Create data file
    filn = f'data.Fe-sc{s}'
    perf = f'data.Fe-sc{s}-perfect'
    os.system(f'/home/c2020spring/share/lammps-data.exe ../Fe-lat-param.txt {s} {s} {s} > {perf}')
    # Modify files for each vacancy
    with open(perf, 'r') as r:
        ls = r.readlines()
        # Vacancy - skip last line
        with open(f'{filn}-V1', 'w') as f:
            for i, line in enumerate(ls):
                if i==2:
                    N = int(line.split(' ')[0])
                    f.write(f'{N-1}  atoms\n')
                elif i < len(ls)-1:
                    f.write(line)
        # Di-vacancy - skip last two lines
        with open(f'{filn}-V2', 'w') as f:
            for i, line in enumerate(ls):
                if i==2:
                    N = int(line.split(' ')[0])
                    f.write(f'{N-2}  atoms\n')
                elif i < len(ls)-2:
                    f.write(line)
        # SIA - modify last line and add one line
        with open(f'{filn}-SIA', 'w') as f:
            for i, line in enumerate(ls):
                if i==2:
                    N = int(line.split(' ')[0])
                    f.write(f'{N+1}  atoms\n')
                elif i < len(ls)-1:
                    f.write(line)
                elif i == len(ls)-1:
                    coord = float(line.rstrip().split(' ')[-1])
                    delta = 0.5
                    f.write(f'{N} 1 {coord-delta} {coord-delta} {coord-delta}\n')
                    f.write(f'{N+1} 1 {coord+delta} {coord+delta} {coord+delta}\n')
    # Create three modified input scripts, launch them and extract the relaxed system energy. 
    ts = ['perfect', 'V1', 'V2', 'SIA']
    for t in ts:
        with open('../input.Fe-opt-defect', 'r') as f:
            inp = f.readlines()
            with open(f'input.Fe-opt-defect-{t}', 'w') as w:
                for line in inp:
                    if 'variable dfile string' in line:
                        w.write(f'variable dfile string data.Fe-sc{s}-{t}\n')
                    else:
                        w.write(line)
        # Launch calculation
        os.system(f'mpirun -np 2 /home/bin/lmp_mpich < input.Fe-opt-defect-{t} > output.Fe-opt-defect-{t}')
        # Extract energy
        with open(f'e.fe', 'r') as ef:
            e = float(ef.readlines()[-2].rstrip().split('=')[1])
            if t == 'perfect':
                Ns = N
            elif t == 'V1':
                Ns = N-1
            elif t == 'V2':
                Ns = N-2
            elif t == 'SIA':
                Ns = N+1
            Ed[f'{s}'][f'{t}'] = {
                                    'E' : e,
                                    'N' : Ns    
                                }
    # Change back dir
    os.chdir('..')
# Print Ed to byte file
with open('Ed.pckl', 'wb') as f:
    pickle.dump(Ed, f)
