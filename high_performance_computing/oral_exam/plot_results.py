# Internal imports
import os

# External imports
import numpy as np
import matplotlib.pyplot as plt

# Set plot params
plt.rc('font', size=18)          # controls default text sizes
plt.rc('axes', titlesize=18)     # fontsize of the axes title
plt.rc('axes', labelsize=18)     # fontsize of the x and y labels
plt.rc('xtick', labelsize=16)    # fontsize of the tick labels
plt.rc('ytick', labelsize=16)    # fontsize of the tick labels
plt.rc('legend', fontsize=14)    # legend fontsize

sizes_original = [1,2,5,10,20,25]
sizes_warp = [1,2,4,8,16,32]
original = ['100x100', '10000x10000', '100000x100']
warp = ['128x128', '10016x10016', '100000x128']

overhead = {}
naive = {}
handin = {}
optimized = {}

for file in os.listdir('bench'):
    fname = file.split('.')[0]
    fsplt = fname.split('_')
    gs1 = fsplt[-2]
    gs2 = fsplt[-1]
    filename = ''
    if len(fsplt) == 3 and not fsplt[0] == 'optimized':
        filename = fsplt[0]
    elif len(fsplt) == 4:
        filename = f'{fsplt[0]}_{fsplt[1]}'
    grid_size = f'{gs1}x{gs2}'
    
    if filename == 'overhead':
        # Open file and read mean and std
        with open(f'bench/{file}', 'r') as f:
            line = f.readlines()[1].split(',')
        mean = float(line[1])
        std = float(line[2])
        overhead[grid_size] = {
            'mean': mean,
            'std': std
        }
    elif filename == 'naive':
        # Open file and read mean and std
        with open(f'bench/{file}', 'r') as f:
            line = f.readlines()[1].split(',')
        mean = float(line[1])
        std = float(line[2])
        naive[grid_size] = {
            'mean': mean,
            'std': std
        }
    elif filename == 'handin':
        # Open file and read mean and std
        with open(f'bench/{file}', 'r') as f:
            line = f.readlines()[1].split(',')
        mean = float(line[1])
        std = float(line[2])
        handin[grid_size] = {
            'mean': mean,
            'std': std
        }
    elif filename == 'optimized_ls':
        # Open file and read mean and std
        local = []
        means = []
        stds = []
        with open(f'bench/{file}', 'r') as f:
            lines = f.readlines()
            for i, line in enumerate(lines[1:]):
                if grid_size in original:
                    ls = line.split(',')
                    if int(ls[-1]) in sizes_original:
                        means.append(float(ls[1]))
                        stds.append(float(ls[2]))
                        local.append(int(ls[-1]))
                elif grid_size in warp:
                    ls = line.split(',')
                    if int(ls[-1]) in sizes_warp:
                        ls = line.split(',')
                        means.append(float(ls[1]))
                        stds.append(float(ls[2]))
                        local.append(int(ls[-1]))
        assert len(means) == len(sizes_original) or len(means) == len(sizes_warp)
        optimized[grid_size] = {
            'mean': means,
            'std' : std,
            'local': local
        }

# Generate plots
porder =  {
    '100x100': [0,0],
    '10000x10000': [0,1],
    '100000x100': [0,2],
    '128x128': [1,0],
    '10016x10016': [1,1],
    '100000x128': [1,2],
    }
fig, axs = plt.subplots(2, 3, figsize=(16,12))
for i, gs in enumerate(optimized):
    opt = optimized[gs]
    p = porder[gs]
    if gs == '100x100' or gs == '10000x10000' or gs == '100000x100':
        ov = overhead[gs]
        na = naive[gs]
        ha = handin[gs]
        axs[p[0]][p[1]].axhline(ov['mean'], c='k', linewidth=2, linestyle='--', label='Overhead')
        axs[p[0]][p[1]].axhline(na['mean'], c='r', linewidth=2, linestyle=':',label='Naive')
        axs[p[0]][p[1]].axhline(ha['mean'], c='g', linewidth=2, linestyle='-.', label='Handin')
    axs[p[0]][p[1]].errorbar(x=opt['local'], y=opt['mean'], yerr=opt['std'], ecolor='b', elinewidth=2, capsize=2, capthick=2, linewidth=2, label='Optimized')
    axs[p[0]][p[1]].set_xlabel('Local size, (n x n)')
    axs[p[0]][p[1]].set_ylabel('Program execution time, s')
    axs[p[0]][p[1]].set_title(f'{gs}')
    axs[p[0]][p[1]].grid()
    axs[p[0]][p[1]].legend(loc='best')
plt.tight_layout()
plt.savefig('local_size.png')

