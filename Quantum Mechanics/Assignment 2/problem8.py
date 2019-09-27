import numpy as np
import matplotlib.pyplot as plt
import random
import mcint


def diff_r_rprim(r, r_p, theta, theta_p, phi, phi_p):
    '''
    Returns the difference between two arbitrary vectors r and r',
    given their position.
    '''
    term1 = r**2 + r_p**2
    term2  = 2*r*r_p*(np.sin(theta)*np.sin(theta_p)*np.cos(theta-phi)+np.cos(theta)*np.cos(theta_p))
    diff_sqrt = np.sqrt(term1-term2)
    return diff_sqrt


def k_dot_r(k,r,theta):
    '''
    Returns the value of the dot product between k and r, 
    with the coordinates aligned with the z || z' axes.
    '''
    return k*r*np.cos(theta)


""" def kp_dot_rp(alpha, beta, kp, rp, thetap, phip):
    '''
    alpha and beta are the azimuthal and the polar angles for
    k' respectively. The integral should be independent of 
    beta, due to spherical symmetry.
    '''
    return kp*rp*(np.sin(alpha)*np.cos(beta-phip)+np.cos(alpha)*np.cos(thetap))


def integrand(x, alpha, beta, k, kp):
    # Dimensions to integrate over
    r = x[0]
    phi = x[1]
    theta = x[2]
    rp = x[3]
    phip = x[4]
    thetap = x[5]
    
    factor1 = np.exp(-1j*kp_dot_rp(alpha, beta, kp, rp, thetap, phip))
    factor2 = np.exp(1j*k_dot_r(k, r, theta))
    diffrrp = diff_r_rprim(r, rp, theta, thetap, phi, phip)
    factor3 = np.exp(1j*k*diffrrp)/diffrrp
    return factor1*factor2*factor3 """

# Kevins code
def differens(r, r_p, theta, theta_p, phi, phi_p):    
    return np.sqrt(r**2 + r_p**2 - 2*r*r_p*(np.sin(theta)*np.sin(theta_p)*np.cos(phi-phi_p) + np.cos(theta)*np.cos(theta_p)))

def dotprod1(k, r_p, alpha, theta_p, beta, phi_p):
    return k*r_p*(np.sin(alpha)*np.sin(theta_p)*np.cos(beta-phi_p) + np.cos(alpha)*np.cos(theta_p))
    
def dotprod2(k,r,theta):
    return k*r*np.sin(theta)
def integrand(k, r, r_p, theta, theta_p, alpha, beta, phi, phi_p):
    exp1 = np.exp(1j*dotprod1(k,r_p,alpha,theta_p,beta,phi_p))
    exp2 = np.exp(1j*dotprod2(k,r,theta))
    exp3 = np.exp(1j*differens(r, r_p, theta, theta_p, phi, phi_p))
    fac = r**2*r_p**2*np.sin(theta)*np.sin(theta_p)/differens(r, r_p, theta, theta_p, phi, phi_p)
    return exp1*exp2*exp3*fac


# constants
alpha = np.linspace(0, 2*np.pi, 20)
beta = 0  # Should be independent of beta
k = 1
kp = k  # Elastic 
R = 1e-10  # Integration limit
n_iter = 10000


# Define sampler
#random.seed(1)  # Remove to check if stable TODO
""" def random_sampler(R):
    r = random.uniform(0., R)
    phi = random.uniform(0, 2.*np.pi)
    theta = random.uniform(0., np.pi)
    rp = random.uniform(0., R)
    phip = random.uniform(0, 2*np.pi)
    thetap = random.uniform(0., np.pi)
    return [r, phi, theta, rp, phip, thetap]


def monte_carlo(n_iter, alpha, beta, k, kp, R):
    # Get random values
    rand = random_sampler(R)
    integ = 0
    for i in range(n_iter):
        # Perform integrand calculation
        integ += integrand(rand, alpha, beta, k, kp)
    return integ/n_iter


# start Monte Carlo integration
sig_alpha = np.zeros(len(alpha))
for idx, a in enumerate(alpha):
    f_alpha = monte_carlo(n_iter, a, beta, k, kp, R)
    sig_alpha[idx] = np.abs(f_alpha)**2 """


N = 10000000
R = 1
k = 1


r     = R*np.random.rand(N)
r_p   = R*np.random.rand(N)
theta = np.pi*np.random.rand(N)
theta_p = np.pi*np.random.rand(N)
phi   = 2*np.pi*np.random.rand(N)
phi_p   = 2*np.pi*np.random.rand(N)

alpha = np.linspace(0,2*np.pi,20)
beta = 0

f_alpha = np.zeros(len(alpha))

for j,alph in enumerate(alpha):
    print(f'{j+1} of {len(alpha)}')
    I =  integrand(k, r, r_p, theta, theta_p, alph, beta, phi, phi_p)
    f_alpha[j] = np.abs(((I)).sum()/N)**2


fig, ax = plt.subplots()
ax.plot(alpha, f_alpha, 'b')
plt.show()