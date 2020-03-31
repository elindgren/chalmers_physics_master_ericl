import scipy.special as spec
import numpy as np
import matplotlib.pyplot as plt


def calc_E_p4(alpha, V_0, zeros):
    E = -zeros/alpha - V_0*np.ones((len(zeros)))
    return E


def calc_E_p3(V_0, hbar, L, m, n=20):
    N = [i+1 for i in range(n)]
    N = np.array(N)
    E = (V_0*hbar/(L*np.sqrt(2*m))*np.pi*(N - 1/4*np.ones((len(N)))*3/2))**(2/3)
    return E


def calc_alpha(m, V_0, hbar, L):
    return (2*m*(L/(V_0*hbar))**2)**(1/3)


# Calculate the 20 first zeros to Airy function
N = 20
ai_z, _, _, _ = spec.ai_zeros(N)

# Constants:
hbar = 1e-34
m = 1e-31
L = 1e-15
V_0 = hbar**2*np.pi**2/(m*L**2)*56169/128
alpha = calc_alpha(m, V_0, hbar, L)


E3 = calc_E_p3(V_0, hbar, L, m, N)
E4 = calc_E_p4(alpha, V_0, ai_z)
fig, ax = plt.subplots()
ax.plot(E3, 'b.')
ax.plot(E4, 'r.')
print("Relative errors:")
print(np.abs(np.divide(E4, E3)))

plt.show()
