import numpy as np
import matplotlib.pyplot as plt


''' 
    Plot the Buckingham potential for the MgO example from the lecture.
    Doesn´t work, but you get the idea. :=)
'''

def coulumb(r12, q1, q2):
    return 14.4*q1*q2/r12

def shortrange(r12, rho, A):
    return A*np.exp(-r12/rho)

def longrange(r12, C):
    return -C/r12**6


def buckingham(r12, q, A, rho, C):
    q1 = q[0]
    q2 = q[1]
    return coulumb( r12, q1, q2 ) + shortrange(r12, rho, A) + longrange(r12, C) 

txt = ['Mg-Mg', 'Mg-O', 'O-O']
A = [0.0, 821.6, 0.342]  # eV
rho = [1e-5, 0.3242, 0.1490]  # Å
C = [0.0, 0.0, 27.88]  # eV Å^6
r_cut = 10.0  # Å

# 
r12 = np.linspace(0.1,r_cut,100)
q1 = 2
q2 = -2
q = [[q1, q1], [q1, q2], [q2, q2]]

fig, ax = plt.subplots(figsize=(8,6))
for i, l in enumerate(txt):
    b = buckingham( r12, q[i], A[i], rho[i], C[i] )
    ax.plot(r12, b, label=l)

ax.set_xlabel('r12')
ax.set_ylabel('Potential eV')
ax.grid()
ax.legend(loc='best')
plt.show()


