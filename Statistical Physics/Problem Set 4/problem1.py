import numpy as np
import matplotlib.pyplot as plt

# Define function

def C(l, t, A, L, R, c):
    k = 2*np.pi*l/L
    return l*(l+1)*A/k**2 * np.abs( np.cos( k * c * np.sqrt(1/(3*(1+R))) * t ) )**2 

# Constants
c = 1  # ly/y
t = 380000  # 380000 years
A = 1
# L = 9.7  # ly
R = 0.7
L = 440 * c * t / np.sqrt((3*(1+R)))
print(f'L={L:.4f} ly')
ls = np.linspace(40, 1500, 1000)

# plot
fig, ax = plt.subplots(figsize=(8,6))
ax.plot(ls, C(ls, t, A, L, R, c))
plt.show()