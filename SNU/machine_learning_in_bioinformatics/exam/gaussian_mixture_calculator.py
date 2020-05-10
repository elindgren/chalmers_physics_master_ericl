# External impors
import numpy as np
import matplotlib.pyplot as plt
from pandas import *

def Nd(x, mu, var):
    sigma = np.sqrt(var)
    return 1/(sigma*np.sqrt(2*np.pi)) * np.exp( -1/2*((x-mu)/sigma)**2 )


def log_like(x, mu, var, pi):
    log_like = 0
    for xn in x:
        log_like += np.log( pi[0]*Nd(xn, mu[0], var[0]) + pi[1]*Nd(xn, mu[1], var[1]) )
    return log_like


def tau_Znk(x, mu, var, pi):
    tZnk = np.zeros((len(x), len(pi) ))
    for n,xn in enumerate(x):
        for k, pik in enumerate(pi):
            ik = 0**k # The other index; when k = 0, ik=1, when k = 1, ik=0
            tZnk[n,k] = pik*Nd(xn, mu[k], var[k])/(pik*Nd(xn, mu[k], var[k]) + pi[ik]*Nd(xn, mu[ik], var[ik]))
    return tZnk


def print_iteration(it, tZnk, Nk, mu, var, pi, logL):
    print(f'******* Iteration {it} *******')
    print('---- E-step ----')
    print(f'tau(Znk):')
    dtZnk = DataFrame(tZnk)
    dtZnk.columns = ['1. M', '2. B']
    dtZnk.index += 1
    print(f'{dtZnk}')
    print()
    print('---- M-step ----')
    print(f'Nk:')
    dNk = DataFrame(Nk)
    dNk.index += 1
    print(f'{dNk}')
    print(f'mu:')
    dmu = DataFrame(mu)
    dmu.index += 1
    print(f'{dmu}')
    print(f'Var:')
    dvar = DataFrame(var)
    dvar.index += 1
    print(f'{dvar}')
    print(f'pi:')
    dpi = DataFrame(pi)
    dpi.index += 1
    print(f'{dpi}')
    print()
    print(f'logL = {logL:.4f}')


# Controls
solve = False  # If true, iterates until convergence
print_step = False
tol = 1e-6

# Data
x = np.array([2,4,7])
N = len(x)
# Initialization
print('******* Initialization *******')
mu = [3, 6]
var = [0.5,0.5]
pi = [0.5, 0.5] 
logL = log_like(x, mu, var, pi)
print(f'logL = {logL:.4f}')
# print(
#     np.log(0.5*Nd(2,3,0.5)+0.5*Nd(2,6,0.5))+
#     np.log(0.5*Nd(4,3,0.5)+0.5*Nd(4,6,0.5))+
#     np.log(0.5*Nd(7,3,0.5)+0.5*Nd(7,6,0.5))
# )
# E-step
tZnk = tau_Znk(x, mu, var, pi)
# print(
#     0.5*Nd(2,3,0.5)/(0.5*Nd(2,3,0.5)+0.5*Nd(2,6,0.5)) 
# ) #z11
# print(
#     0.5*Nd(4,6,0.5)/(0.5*Nd(4,3,0.5)+0.5*Nd(4,6,0.5))
# ) #z22

# M-step 
Nk = np.sum(tZnk, axis=0)
# print(
#     9.999997e-01+9.525741e-01+3.059022e-07
# )
mu = np.sum(tZnk.T*x,axis=1)/Nk
# print(
#     (2*9.999997e-01+4*9.525741e-01+7*3.059022e-07)/1.952574
# )
# print(
#     (2*3.059022e-07+4*4.742587e-02+7*9.999997e-01)/1.047426
# )
var = np.sum(np.array( [tZnk[n,:]*(xn-mu)**2 for n, xn in enumerate(x)] ), axis=0)/Nk
# print(
#     (9.999997e-01*(2-2.975712)**2+9.525741e-01*(4-2.975712)**2+3.059022e-07*(7-2.975712)**2)/1.952574
# )
# print(
#     (3.059022e-07*(2-6.864163)**2+4.742587e-02*(4-6.864163)**2+9.999997e-01*(7-6.864163)**2)/1.047426
# )
pi = Nk/N
# print(
#     1.952574/3
# )
# print(
#     1.047426/3
# )

# Calculate likelihood
logL_old = logL
logL = log_like(x, mu, var, pi)
# print(
#     np.log(0.650858*Nd(2,2.975712,0.999412)+0.349142*Nd(2,6.864163,0.389062))+
#     np.log(0.650858*Nd(4,2.975712,0.999412)+0.349142*Nd(4,6.864163,0.389062))+
#     np.log(0.650858*Nd(7,2.975712,0.999412)+0.349142*Nd(7,6.864163,0.389062))
# )

print_iteration(1, tZnk, Nk, mu, var, pi, logL)


# Iterate until convergence

if(solve):
    it = 2
    logLs = [logL_old, logL]
    while(np.abs(logL-logL_old) > tol):
        # E-step
        tZnk = tau_Znk(x, mu, var, pi)
        # M-step
        Nk = np.sum(tZnk, axis=0)
        mu = np.sum(tZnk.T*x,axis=1)/Nk
        var = np.sum(np.array( [tZnk[n,:]*(xn-mu)**2 for n, xn in enumerate(x)] ), axis=0)/Nk
        pi = Nk/N
        # Likelihood
        logL_old = logL
        logL = log_like(x, mu, var, pi)
        logLs.append(logL)
        if print_step:
            print_iteration(it, tZnk, Nk, mu, var, pi, logL)
        if(it > 100):
            break
        it += 1    
        
    # Plot convergence
    fig, ax = plt.subplots(figsize=(8,6))
    ax.plot(logLs)
    ax.set_xlabel('Iterations')
    ax.set_ylabel(r'log likelihood $\log(p(x | \mu, \Sigma, \pi))$')
    plt.show()



