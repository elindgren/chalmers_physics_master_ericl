import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

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