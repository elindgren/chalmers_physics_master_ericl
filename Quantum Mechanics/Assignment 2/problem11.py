import numpy as np
import matplotlib.pyplot as plt

# b)
def a_over_b(gammab):
    fun = -(np.tan(gammab)**(-1)-1-1/gammab)
    return fun


t = np.linspace(0.0000001,10,10000)

fig, ax = plt.subplots()
ax.plot(t, a_over_b(t))
ax.set_ylim(-25,25)
plt.grid()
plt.show()

# c)
# Ramsauer-Townsend effect. tan(delta0)/k = 0 => delta0 = n*pi =>
# => f(theta) = 0 for these energies

# find smallest value of gammab TODO

# d) infinites for gammab = n*pi (n=0,1,2) and gammab=0

# e)