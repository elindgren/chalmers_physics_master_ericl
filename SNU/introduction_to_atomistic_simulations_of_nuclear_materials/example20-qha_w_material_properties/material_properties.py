# External imports
import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_table('w.output', sep='    ', skiprows=2, header=None, names=['T', 'E_0', 'B_0', "B'_0", 'V0'])


# Plot bulk modulus as a function of temperature
fig, ax = plt.subplots(figsize=(8,6))
ax.plot(df['T'], df['B_0'], linewidth=2)
ax.grid()
ax.set_xlabel('Temperature, (K)')
ax.set_ylabel('Bulk modulus, (GPa)')
plt.show()