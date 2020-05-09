# External imports
import numpy as np
import matplotlib.pyplot as plt
from pandas import *    


def split_sequences(S, W):
    X = []
    for s in S:
        X.append( [s[i:i+W] for i in range(len(s)-W+1)] )
    X = np.array( X )
    X = X.flatten()
    return X


def translate_alphabet(x):
    # Translates a sequence from alphabet (string) to an integer array
    x_ser = []
    for ai in x:
        if ai=='A':
            x_ser.append(0)
        elif ai=='C':
            x_ser.append(1)
        elif ai=='G':
            x_ser.append(2)
        elif ai=='T':
            x_ser.append(3)
    x_ser = np.array( x_ser )
    return x_ser


def count_occurences(x_seq, W):
    # Count the number of counts for the different alphabet characters where x_seq is a translated string
    counts = []
    for a in range(W):
        counts.append(len( np.where(x_seq==a)[0] ))
    counts = np.array( counts )
    return counts


def initialize_f(X, W, L, motif=0):
    f = np.zeros((W+1, L))
    n = len(X)
    p = 0.4 # Initialization probability if a character is found at index p in motif
    # Chose the first sequence to be motif
    mot = X[motif]
    bgd = np.delete(X, motif)

    # Background counts for fij
    for b in bgd:
        b = translate_alphabet(b)
        counts = count_occurences(b, W)
        f[0,:] += counts    
    f[0,:] /= (n-1)*W
    assert np.sum(f[0,:]) == 1.0

    # Motif parameter estimation
    mot = translate_alphabet(mot)
    for i, ai in enumerate(mot):
        for j in range(L):
            if  ai==j:
                f[i+1,j] = p
            else:
                f[i+1,j] = (1-p)/(W-1)
    assert np.all(np.sum(f[1:, :], axis=1)-1.0 < 1e-10)
    return f


def model(x, f, bg_model):
    # Calculates p(x|theta) as a multinomial probability
    x_seq = translate_alphabet(x)
    if bg_model:
        p = np.prod(np.array( [f[0, ai] for ai in x_seq] ))
    else: 
        p = np.prod(np.array( [f[i+1, ai] for i, ai in enumerate(x_seq)] ))
    return p


def calc_Z(X, f, lam):
    Z = np.zeros((len(X), len(lam)))
    for i, xi in enumerate(X):
        for j, lamj in enumerate(lam):
            Z[i,j] = model(xi, f, j)*lamj / (model(xi, f, 0)*lam[0] + model(xi, f, 1)*lam[1])
    return Z


def log_likelihood(X, Z, f, lam):
    logL = 0
    for i in range(Z.shape[0]):
        for j in range(Z.shape[1]):
            logL += Z[i,j]*np.log(model(X[i], f, bg_model=j)*lam[j])
    return logL


def find_motif(X, f):
    probs = []
    for x in X:
        probs.append( model(x, f, bg_model=0) )  # prob of sequence x to be motif
    return X[np.argmax(probs)], np.max(probs)


def calc_count(X, Z, W, L):
    c = np.zeros((W+1, L))
    # Background model
    for k in range(L):
        for i, xi in enumerate(X):
            xi = translate_alphabet(xi)
            for j, xij in enumerate(xi):
                if xij==k:
                    c[0, k] += Z[i,1]
    
    # Motif model
    for k in range(L):
        for i, xi in enumerate(X):
            xi = translate_alphabet(xi)
            for j in range(1,W+1):
                xij = xi[j-1]
                if xij==k:
                    c[j, k] += Z[i,0]
    return c


def update_f(c, b):
    f = np.zeros(c.shape)
    for i in range(c.shape[0]):
        for j in range(c.shape[1]):
            f[i,j] = (c[i,j] + b[j]) / (c[i,:].sum() + b.sum())
    return f


# Controls
solve = True  # If true, iterates until convergence
tol = 1e-6

# Sequences
S = ['ATTGGC', 'TTGACC']

# Initialization
lam = [0.5, 0.5]  # Lambda
W = 4  # Motif length
L = 4  # Alphabet length
b = np.ones(W)  # Pseudocounts
X = split_sequences(S, W)
n = len(X)

print('******* Initialization *******')
Xd = DataFrame(X)
Xd.index += 1
print(f'X:')
print(f'{Xd}')
print()
# Initialize f weight matrix
f = initialize_f(X, W, L, motif=0)
print(f'fij:')
print(f'{DataFrame(f)}')
print()

print('******* Algorithm *******')
# E-step: calculate latent variables zij for each sample xi where j is either 0 or 1 (motif or non-motif)
print('---- E-step ----')
Z = calc_Z(X, f, lam)
zd = DataFrame(Z)
zd.index += 1 # shift indices to fit paper convention
zd.columns = [1, 2]
print(f'Zij:')
print(f'{zd}')
print()
# Calculate log-likelihood
logL = log_likelihood(X, Z, f, lam)
print(f'log likelihood: {logL:.8f}')
print()

# M-step
print('---- M-step ----')
c = calc_count(X, Z, W, L)
print(f'cij:')
print(f'{DataFrame(c)}')
print()

f = update_f(c, b)
print(f'updated f:')
print(f'{DataFrame(f)}')
print()

lam = np.sum( Z, axis=0 ) / n
print(f'updated lambda:')
print(f'{DataFrame(lam)}')
print()

# Calculate log-likelihood
logL_old = logL
logL = log_likelihood(X, Z, f, lam)
print(f'log likelihood: {logL:.8f}')
print()


##### Converge

# Iterate until convergence

if(solve):
    it = 1
    logLs = [logL_old, logL]
    while(np.abs(logL-logL_old) > tol):
        # E-step
        Z = calc_Z(X, f, lam)
        # M-step
        c = calc_count(X, Z, W, L)
        f = update_f(c, b)
        lam = np.sum( Z, axis=0 ) / n
        # Likelihood
        logL_old = logL
        logL = log_likelihood(X, Z, f, lam)
        logLs.append(logL)
        it += 1
        if(it > 100):
            break
    
    # Plot convergence
    fig, ax = plt.subplots(figsize=(8,6))
    ax.plot(logLs)
    ax.set_xlabel('Iterations')
    ax.set_ylabel(r'log likelihood $\log(p(x | \mu, \Sigma, \pi))$')
    # plt.show()


# Find the most probable sequence - our motif
motif, p = find_motif(X, f)
print(f'motif: {motif}, p={p:.4f}')
print()