import numpy as np
import matplotlib.pyplot as plt


# Define constants
hbar = 1
m = 1
# Parameters
V0 = 1
a = 1
L = 1

# Define energy space
E = np.linspace(0, 10, 100)

# Define function
fun = np.sqrt(2*m)/hbar * (-L/V0) * 2/3 * (E**(2/3) - (E+V0)**(2/3))
# Target points: pi/4 + n*pi
n = np.array(range(20))  # 20 target points

target = np.pi/4*np.ones((len(n),)) + np.pi*n
# Plot function

fig, ax = plt.subplots()
#ax.vline(target, [0,1], 'b')
