# External impors
import numpy as np
import matplotlib.pyplot as plt
from pandas import *

def Nd(x, mu, sigma):
    return 1/(sigma*np.sqrt(2*np.pi)) * np.exp( -1/2*((x-mu)/sigma)**2 )


def log_like(x, mu, sigma, pi):
    log_like = 0
    for xn in x:
        log_like += np.log( pi[0]*Nd(xn, mu[0], sigma[0]) + pi[1]*Nd(xn, mu[1], sigma[1]) )
    return log_like


def tau_Znk(x, mu, sigma, pi):
    tZnk = np.zeros((len(x), len(pi) ))
    for n,xn in enumerate(x):
        for k, pik in enumerate(pi):
            ik = 0**k # The other index; when k = 0, ik=1, when k = 1, ik=0
            tZnk[n,k] = pik*Nd(xn, mu[k], sigma[k])/(pik*Nd(xn, mu[k], sigma[k]) + pi[ik]*Nd(xn, mu[ik], sigma[ik]))
    return tZnk


# Controls
solve = False  # If true, iterates until convergence
tol = 1e-6

# Data
x = np.array([1, 5, 9, 11.5, 3])
N = len(x)
# Initialization
print('******* Initialization *******')
mu = [2, 5]
sigma = [1,1]
pi = [0.25, 0.75]
logL = log_like(x, mu, sigma, pi)
print(f'logL = {logL:.4f}')


print('******* Iteration 1 *******')
# E-step
print('---- E-step ----')
tZnk = tau_Znk(x, mu, sigma, pi)
print(f'tau(Znk):')
print(f'{DataFrame(tZnk)}')
print()
# M-step 
print('---- M-step ----')
Nk = np.sum(tZnk, axis=0)
print(f'Nk:')
print(f'{DataFrame(Nk)}')

mu = np.sum(tZnk.T*x,axis=1)/Nk
print(f'mu:')
print(f'{DataFrame(mu)}')

sigma = np.sum(np.array( [tZnk[n,:]*(xn-mu)**2 for n, xn in enumerate(x)] ), axis=0)/Nk
print(f'Sigma:')
print(f'{DataFrame(sigma)}')
# print((9.9834e-1*(1-1.7575)**2 + 3.6893e-3*(5-1.7575)**2 + 2.2752e-8*(9-1.7575)**2 + 1.2584e-11*(11-1.7575)**2 + 5.9902e-1*(3-1.7575)**2)/1.6011)

pi = Nk/N
print(f'pi:')
print(f'{DataFrame(pi)}')
print()

# Calculate likelihood
logL_old = logL
logL = log_like(x, mu, sigma, pi)
print(f'logL = {logL:.4f}')

# Iterate until convergence

if(solve):
    it = 1
    logLs = [logL_old, logL]
    while(np.abs(logL-logL_old) > tol):
        # E-step
        tZnk = tau_Znk(x, mu, sigma, pi)
        # M-step
        Nk = np.sum(tZnk, axis=0)
        mu = np.sum(tZnk.T*x,axis=1)/Nk
        sigma = np.sum(np.array( [tZnk[n,:]*(xn-mu)**2 for n, xn in enumerate(x)] ), axis=0)/Nk
        pi = Nk/N
        # Likelihood
        logL_old = logL
        logL = log_like(x, mu, sigma, pi)
        logLs.append(logL)
        it += 1
        if(it > 100):
            break
    
    # Plot convergence
    fig, ax = plt.subplots(figsize=(8,6))
    ax.plot(logLs)
    ax.set_xlabel('Iterations')
    ax.set_ylabel(r'log likelihood $\log(p(x | \mu, \Sigma, \pi))$')
    plt.show()



