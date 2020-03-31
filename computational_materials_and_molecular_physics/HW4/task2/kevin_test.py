import numpy as np
import matplotlib.pyplot as plt


### from helper functions
def fold(x_t, x_i, y_i, width):
    '''
        Convolutes each peak in the discrete spectrum
        with a Gaussian of chosen width.
        
        inputs:
        x_t:    vector with energies (i.e. linspace/np.arange)
        x_i:    stick spectrum energies
        y_i:    stick spectrum intensities
        width:  width of the Gaussian

        outputs:
        y_t:    convoluted spectrum, i.e. intensities
    '''
    def Gauss(x0):
        norm = 1.0 / (width * np.sqrt(2 * np.pi))
        return norm * np.exp(-0.5 * (x_t - x0)**2 / width**2)

    y_t = np.zeros_like(x_t)
    for x, y in zip(x_i, y_i):
        y_t += y * Gauss(x)
    return y_t


def calc_omega(a = 1):
    dump = np.load('./dumpTask1.npz')
    i_p = dump['i_p']
    a_p = dump['a_p']
    w_p = dump['ediff_p']
    n_p = dump['fdiff_p']
    mu_x = dump['mux_p']
    mu_y = dump['muy_p']
    mu_z = dump['muz_p']
    K = dump['K_pp']

    N = len(i_p)
    # omega = np.zeros((N,N))
    
    test = np.sqrt(w_p*n_p)
    omega = np.diag(w_p**2, k=0) + 2*np.outer(test,test)*K*a

    # Omega = np.diag(v=w_p**2, k=0) # First add diagonal
    # for p in range(len(Omega)):
    #     for q in range(len(Omega)):
    #         Omega[p,q] += 2*np.sqrt(n_p[p]*w_p[p]) * K[p,q] * np.sqrt(n_p[q]*w_p[q])

    #alternativt sätt samma resultat
    #for p in range(N):
    #    for q in range(N):
    #        omega[p,q] = 2*np.sqrt(n_p[p]*w_p[p]*w_p[q]*n_p[q])*K[p,q]*a
    #        if p == q:
    #           omega[p,q] += w_p[p]**2

    return omega,n_p

def calc_photo_abs(omega):
    dump = np.load('./dumpTask1.npz')
    i_p = dump['i_p']
    a_p = dump['a_p']
    w_p = dump['ediff_p']
    n_p = dump['fdiff_p']
    mu_x = dump['mux_p']
    mu_y = dump['muy_p']
    mu_z = dump['muz_p']
    K = dump['K_pp']

    
    eig_val, eig_vec = np.linalg.eig(omega)
    
    f = np.zeros((len(eig_val),3))
    mu = np.zeros((len(eig_val),3))
    mu[:,0] = mu_x
    mu[:,1] = mu_y
    mu[:,2] = mu_z

    # eigVal, F = np.linalg.eig(omega)  # Eigenvalue, eigenvectors
    # sort = np.argsort(eigVal)
    # eigVal = eigVal[sort]
    # F = F[:,sort]

    # # Extract the excitations in Ha, and convert them to eV
    sortOrder = np.argsort(eig_val)
    eig_val = eig_val[sortOrder]
    eig_vec = eig_vec[:,sortOrder]
    eps = np.sqrt(eig_val[sortOrder]) * 27.2

    # # Calculate the oscillator strength
    # f = np.zeros(len(eigVal))
    # for i, e in enumerate(eigVal):
    #     for alpha, mua_p in enumerate(mu_p):
    #         f_ia = 2 * np.abs( np.sum( mua_p * np.sqrt(n_p*w_p) * F[:,i] ) )
    #         f[i] += f_ia
    #     f[i] /= 3 # Average over all dipole moments


    for x in range(3):
        for I in range(len(eig_val)):
            f[:,I] += np.sum(eig_vec[:,I]* mu[p,x] *np.sqrt(w_p[p]*n_p[p]))
            
    f = 2*(np.abs(f))**2
    print(f.shape)

    return f, eps

lwidth=2
fsize=14

omega,n_p = calc_omega()
f, energy = calc_photo_abs(omega)

### plot descrite spectrum 
fig,ax = plt.subplots(figsize = (8,4))
ax.plot(energy,f.mean(axis = 1))
#ax.plot(energy,f)

### convolut and plot spectrum 
w = np.linspace(1,6.2,10000)
a = fold(w, energy, f.mean(axis = 1), 0.06)

fig,ax = plt.subplots(figsize = (8,4))
ax.plot(w, a, Linewidth = lwidth, label = "Casida problem 2")
# ax.plot(spec_p1[:,0],spec_p1[:,1],"--" ,Linewidth = lwidth, label = "Casida problem 1")

ax.legend(fontsize = fsize)
ax.set_xlabel("$\omega$ [eV]", fontsize = fsize) # This units is not correct ? 
ax.set_ylabel(r"$Photoabsorption$ [1/eV]", fontsize = fsize)   # This units is not correct ? 
ax.tick_params(labelsize=fsize)
plt.show()