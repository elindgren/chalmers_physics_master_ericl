import numpy as np
import matplotlib.pyplot as plt

def calc_q(k, alpha):
    return 2*k*np.sin(alpha/2)


def calc_beta(k,R):
    beta = 1/k*(6/k**4+np.exp(-1j*k*R)*(1j*R**3/k - 3*R**2/k**2 - 1j*6*R/k**3 - 6/k**4)) + 2/k**2*(np.sin(k*R)-R**2/k*np.cos(k*R))*(np.exp(1j*k*R)*(-1j*R**2/k + R/k**2 + 2j/k**3)-2j/k**3) - 2/k**2*(2*R/k**2*np.sin(k*R)-R**2/k*np.cos(k*R)+2/k**3*np.cos(k*R) + 2/k**3) 
    #beta = 1/k*(6/k**4+np.exp(-1j*k*R))
    return beta

def calc_f(alpha, k, R, m, v0, hbar):
    q = calc_q(k, alpha)
    term1 = 1/q**3 * (np.sin(q*R)-R*q*np.cos(q*R))
    term2 = 1/q**4 * ((2-R**2*q**2)*np.cos(q*R)+2*q*R*np.sin(q*R)-2)
    factor = -1/(4*np.pi)*2*m/hbar**2
    return factor*(term1-term2)



# constants
hbar = 4.136e-15  # ev s
R = 1e-10  # 1 ångström
v0 = 1  # 1 ev
E = 10 # 1 ev
m = 938e6  # 938 MeV/c^2

k = np.sqrt(2*m*E)/hbar
eps = 1e-2
alpha = np.linspace(0+eps,2*np.pi-eps,300)

dsig_prob6 = [(calc_f(a, k, R, m, v0, hbar))**2 for a in alpha]


fig, ax = plt.subplots()
ax.plot(alpha, dsig_prob6, 'b')
ax.set_ylim(0,1e-24)
ax.set_xlabel("alpha, from 0 to 2*pi")
ax.set_ylabel("differential cross section")
plt.savefig("problem6.png")
plt.title("First order Born approximated diff cross section")

plt.show()