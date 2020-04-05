# External imports
import numpy as np
import matplotlib.pyplot as plt


# Set plot params
plt.rc('font', size=18)          # controls default text sizes
plt.rc('axes', titlesize=18)     # fontsize of the axes title
plt.rc('axes', labelsize=18)     # fontsize of the x and y labels
plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
plt.rc('ytick', labelsize=18)    # fontsize of the tick labels
plt.rc('legend', fontsize=18)    # legend fontsize

'''
    Implement EM for a Gaussian Mixture problem with two clusters. 
'''

def gauss(x, mu, sigma):
    return 1/np.sqrt(2*np.pi*sigma**2) * np.exp( -(x-mu)**2 / (2*sigma**2) )

def gen_data(N=10000):
    '''
        Generate a data set of size N, where each data point 
        has a prob of coming from one of two Gaussians.
    '''
    f = 0.3
    d_g1 = np.random.normal(loc=10, scale=2, size=int( f*N ))
    d_g2 = np.random.normal(loc=20, scale=5, size=int( (1-f)*N ))
    return np.hstack((d_g1, d_g2)), d_g1, d_g2

def em(d, plot=True):
    # Initialization
    # Gauss 1
    mu_1 = np.random.rand()*10
    sigma_1 = np.random.rand()*5
    phi_1 = 0.5 * np.ones(len(d))  # probability of a datapoint i coming from gaussian 1
    # Gauss 2
    mu_2 = np.random.rand()*10
    sigma_2 = np.random.rand()*5
    phi_2 = 1 - phi_1
    if plot:
        fig, ax = plt.subplots(figsize=(8,6))

    mu_1_old = 0
    tol = 1e-4
    i = 0
    # Recurrence
    while np.abs(mu_1 - mu_1_old) > tol:
        mu_1_old = mu_1
        # Estimate latent variables given current parameters
        l_1 = gauss(d, mu_1, sigma_1)  # likelihood vector for data belonging to gauss 1
        l_2 = gauss(d, mu_2, sigma_2) 

        # Get posterior probabilities for the data belonging to each gaussian, including prior phi
        b_1 = l_1*phi_1 / (l_1*phi_1 + l_2*phi_2)
        b_2 = l_2*phi_2 / (l_1*phi_1 + l_2*phi_2)

        # Calculate expected value of new parameters - MLE (expectation value)
        mu_1 = np.sum( b_1*d ) / np.sum( b_1 )
        sigma_1 = np.sqrt( np.sum(b_1*(d-mu_1)**2) / np.sum(b_1) )
        phi_1 = 1/len(d) * np.sum(b_1)

        mu_2 = np.sum( b_2*d ) / np.sum( b_2 )
        sigma_2 = np.sqrt( np.sum(b_2*(d-mu_2)**2) / np.sum(b_2) )
        phi_2 = 1/len(d) * np.sum(b_2)

        # Plot 
        if plot and i%10==0:
            x = np.linspace(0, 40, 500)
            ax.plot(x, gauss(x, mu_1, sigma_1))
            ax.plot(x, gauss(x, mu_2, sigma_2))
        i += 1
    print(f'EM converged in {i} iterations')
    return [mu_1, mu_2], [sigma_1, sigma_2], [phi_1, phi_2]

np.random.seed(1)

d, d_g1, d_g2 = gen_data(10000)

mu, sigma, phi = em(d, plot=False)

# Plot final distributions
fig, ax = plt.subplots(figsize=(8,6))
x = np.linspace(0, 40, 500)
ax.plot(x, gauss(x, mu[0], sigma[0]), linewidth=2, label='Gauss 1')
ax.plot(x, gauss(x, mu[1], sigma[1]), linewidth=2, label='Gauss 2')

# Plot data distributions
ax.hist(d, bins=100, density=True, label='True data')
ax.hist(d_g1, bins=100, density=True, label=r'd_{G1}')
ax.hist(d_g2, bins=100, density=True, label=r'd_{G2}')
ax.legend(loc='best')
ax.grid()
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'p(x)')
ax.set_title(r'EM for GMM, $k=2$')
plt.tight_layout()
plt.savefig('em_gmm_k=2.png')


