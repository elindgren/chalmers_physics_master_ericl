# Internal package 
import pprint

# External packages
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
    Python implementation for the analysis of the Dishonest Casino example from
    Homework 1 in the course Machine Learning in Bioinformatics, Spring 2020. 

    All code is given in this single file for comprehensive purposes.

    Created by Eric Lindgren, SNU ID ericlin
    Spring 2020, Seoul
'''

def viterbi(x, a, e):
    '''
        Compute the most probable state path pi* and the probability of that path p(x,pi*)
        for a given emission sequence, given HMM transition parameters a and emission probabilities b. 

        Perform computations in log-space for the sake of numerical precision
    '''
    v = np.zeros((len(x)+1, 3)) # v matrix in real space - 1 longer than x since start in state B
    ptr = np.zeros((len(x), 3))  # pointer for most probable state
    pi = np.zeros(len(x), dtype=np.int)
    pi_test = np.zeros(len(x), dtype=np.int)

    #**** Initialisation ****
    v[0,:] = [1, 0, 0]  # Start in begin state
    # Convert to log-space
    with np.errstate(divide='ignore'):
        V = np.log( v )
        A = np.log( a )  # Note! Not count variable A, log-space equivalent of a!
        E = np.log( e )
    
    #**** Recursion ****
    for i in range(1, len(x)+1):
        xi = x[i-1]
        for l in range(1, V.shape[1]):
            # Can't go back to B, so skip that case
            VA = V[i-1,:] + A[:,l]
            V[i,l] = E[l-1, xi] + np.max( VA )
            ptr[i-1, l] = np.argmax( VA )  #! Arguments off by 1
        pi_test[i-1] = np.argmax( V[i,:] )  #! Off by one? # For assertion purposes

    #**** Termination ****
    a_k0 = [0, 1, 1]  # a_k0: probability of going from any state to end state is 1!
    with np.errstate(divide='ignore'):
        A_k0 = np.log( a_k0 )
    p_pxpi = np.exp(np.max( V[-1,:] + A_k0 ))  # No a since no transition to end state
    pi[-1] = np.argmax( V[-1,:] + A_k0 )  # Exp is monotonous, no need to transform back

    #**** Traceback ****
    for i in reversed(range(1, len(x))):
        pi[i-1] = ptr[i, pi[i]]
    assert np.array_equal(pi, pi_test)

    print(f'Viterbi algorithm: P(x, pi*) = {p_pxpi:.4e}')
    return pi, p_pxpi


def forward(x, a, e):
    '''
        Compute the probability for the sequence up to and including xi in which the model is in state k, given in the matrix f. 
        Also compute the likelihood of the observation, p(x).

        Perform computations in log-space for the sake of numerical precision
    '''
    f = np.zeros((len(x)+1, 3)) # v matrix in real space - 1 longer than x since start in state B

    #**** Initialisation ****
    f[0,:] = [1, 0, 0]  # Start in begin state
    # Convert to log-space
    with np.errstate(divide='ignore'):
        F = np.log( f )
        A = np.log( a )  # Note! Not count variable A, log-space equivalent of a!
        E = np.log( e )
    
    #**** Recursion ****
    for i in range(1, len(x)+1):
        xi = x[i-1]
        for l in range(1, F.shape[1]):
            # Can't go back to B, so skip that case
            FA = F[i-1,:] + A[:,l]
            F[i,l] = E[l-1, xi] + np.log(np.sum(np.exp( FA )))  #* Summation must be performed in real space!
    #**** Termination ****
    a_k0 = [0, 1, 1]  # a_k0: probability of going from any state to end state is 1!
    p = np.sum(np.exp( F[-1,:] )*a_k0)  # No a since no transition to end state
    f = np.exp( F )
    print(f'Forward algorithm: P(x) = {p:.4e}')
    return f, p


def backward(x, a, e):
    '''
        Compute the probability for the sequence end, given that the model is in state k at step i, given in the matrix b. 
        Also compute the likelihood of the observation, p(x).

        Perform computations in log-space for the sake of numerical precision
    '''
    b = np.zeros((len(x)+1, 3)) # v matrix in real space - 1 longer than x since start in state B

    #**** Initialisation ****
    a_k0 = [0, 1, 1]  # a_k0: probability of going from any state to end state is 1!
    b[-1,:] = a_k0
    # Convert to log-space
    with np.errstate(divide='ignore'):
        B = np.log( b )
        A = np.log( a )  # Note! Not count variable A, log-space equivalent of a!
        E = np.log( e )
    
    #**** Recursion ****
    for i in reversed(range(1, len(x))):
        xi = x[i-1]
        for l in range(1, B.shape[1]):
            # Can't go back to B, so skip that case
            BA = B[i+1,:] + A[:,l]
            B[i,l] = np.log(np.sum(np.exp( BA + E[l-1, xi] )))  # Summation must be performed in real space!
    #**** Termination ****
    p = np.sum(np.exp( B[1,1:] + A[0,1:] + E[:,x[0]] ))  # Skip transition a_00, since it's illegal
    b = np.exp( B )
    print(f'Backward algorithm: P(x) = {p:.4e}')
    return b, p

#**** Given parameters, from homework description ****
x_raw = 'HHHTHHTHHHHTTHHHHHHTHHTHHTHHHTHHTTHHHTHHHHHT'  # Raw observation
a = np.array([
    [0, 0.5, 0.5], 
    [0, 0.8, 0.2],
    [0, 0.3, 0.7]
])  # transition matrix a_kl. k=0 is B, k=1 is F, k=2 is L.

e = np.array([
    [0.5, 0.5],
    [0.9, 0.1],
])  # emission probability matrix e_k(b) for states F=0 and L=1. Column 0 is H, 1 is T 

#****  Pre-processing: replace H with 0 and T with 1 in x ****
x = np.array([ int(s=='H') for s in x_raw ])

#****  Obtain most probable state path: Viterbi algorithm ****
pi_m, p_pxpi = viterbi(x, a, e)

#**** Posterior decoding - obtain f and b matrices ****
f, p_f = forward(x, a, e)
b, p_b = backward(x, a, e)
assert np.abs(p_f-p_b)/p_f < 0.15  # TODO The error is around 15%
p = (p_f + p_b) / 2  # Average probability
p_post = f*b/p  # Posterior probability
# Extract highest probability sequence - not necessarily allowed transitions-wise!
pi_hat = np.argmax(p_post, axis=1)

#**** Compare Viterbi and posterior ****
pp = pprint.PrettyPrinter(indent=4)
print('****** Viterbi: pi_m ******')
pp.pprint(pi_m)  #! Probably wrong
print('****** Posterior: pi_hat ******')
pp.pprint(pi_hat)
print('****** Sum of posterior probabilities Sum_k( p(xi=k|x) ) ******')
pp.pprint(np.sum(p_post, axis=1))  # TODO Check all probabilities sum to 1
print('****** Posterior probaility p(xi=k|x) ******')
pp.pprint(p_post)