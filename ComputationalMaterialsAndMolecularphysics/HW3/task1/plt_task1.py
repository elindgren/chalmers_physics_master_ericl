import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

data = np.loadtxt('log_mdTask1.txt', skiprows=1)

fig, ax = plt.subplots(figsize=(8,6))
# data.plot(kind='line', x='Time', y='T', ax=ax)
ax.plot(data[:,0], data[:,4])
ax.axhline(350)
print(data)
plt.show()