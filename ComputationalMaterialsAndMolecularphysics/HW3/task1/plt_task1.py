import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Set plot params
plt.rc('font', size=14)          # controls default text sizes
plt.rc('axes', titlesize=14)     # fontsize of the axes title
plt.rc('axes', labelsize=16)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=14)    # fontsize of the tick labels
plt.rc('ytick', labelsize=14)    # fontsize of the tick labels
plt.rc('legend', fontsize=14)    # legend fontsize

data = np.loadtxt('log_mdTask1.txt', skiprows=1)
avgT = np.mean(data[:,4])

fig, ax = plt.subplots(figsize=(8,6))
# data.plot(kind='line', x='Time', y='T', ax=ax)
ax.plot(data[:,0], data[:,4], linestyle='-', label=r'$T(t)$')
ax.axhline(avgT, c='k', linewidth=2, label=f'Avg. T = {avgT:.2f} K')
ax.set_xlabel('Time, (ps)')
ax.set_ylabel(r'$T$, (K)')
ax.legend(loc='best')
ax.grid()
plt.savefig('task1_T_traj.png')
plt.show()