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
    ptr = np.zeros((len(x)+1, 3))  # pointer for most probable state
    pi = np.ones(len(x)+1, dtype=np.int)  # Include the begin state in the most probable path

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
            V[i,l] = E[l-1, xi] + np.max( VA )  # E[l-1] since we skip state B
            ptr[i , l] = np.argmax( VA )
    #**** Termination ****
    a_k0 = [0, 1, 1]  # a_k0: probability of going from any state to end state is 1!
    with np.errstate(divide='ignore'):
        A_k0 = np.log( a_k0 )
    p_pxpi = np.exp(np.max( V[-1,:] + A_k0 ))  # No a since no transition to end state
    pi[-1] = np.argmax( V[-1,:] + A_k0 )  # Exp is monotonous, no need to transform back

    #**** Traceback ****
    for i in reversed(range(1, len(x)+1)):
        pi[i-1] = ptr[i, pi[i]]
    pi = pi[1:]  # Only return the part that corresponds to the sequence - skip begin state
    
    v = np.exp( V )
    print(f'Viterbi algorithm: P(x, pi*) = {p_pxpi:.4e}')
    return v, pi, p_pxpi


def forward(x, a, e):
    '''
        Compute the probability for the sequence up to and including xi in which the model is in state k, given in the matrix f. 
        Also compute the likelihood of the observation, p(x).

        Return only f values being connected to x: f[1:] - Skip begin state

        Perform computations in log-space for the sake of numerical precision
    '''
    f = np.zeros((len(x)+1, 3)) # v matrix in real space - 1 longer than x for start state

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
    with np.errstate(divide='ignore'):
        A_k0 = np.log( a_k0 )
    p = np.sum(np.exp( F[-1,:] + A_k0 ))   # No a since no transition to end state
    f = np.exp( F )
    print(f'Forward algorithm: P(x) = {p:.20e}')
    return f[1:], p 


def backward(x, a, e):
    '''
        Compute the probability for the sequence end, given that the model is in state k at step i, given in the matrix b. 
        Also compute the likelihood of the observation, p(x).

        Perform computations in log-space for the sake of numerical precision
    '''
    b = np.zeros((len(x), 3)) # v matrix in real space - 1 longer than x since end state
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
        xi = x[i]
        for k in range(1, B.shape[1]):
            # Can't go back to B, so skip that case
            BA = B[i,1:] + A[k,1:] + E[:, xi]
            B[i-1, k] = np.log(np.sum( np.exp( BA )))  # Summation must be performed in real space!
    #**** Termination ****
    p = np.sum(np.exp( A[0,1:] + B[0,1:] + E[:,x[0]] ))  # Skip transition a_00, since it's illegal
    b = np.exp( B )
    print(f'Backward algorithm: P(x) = {p:.20e}')
    return b, p


def translate_to_state(pi):
    ''' Convert numerical state sequence to F/L '''
    return np.array( ['F' if s==1 else 'L' for s in pi] )


#**** Given parameters, from homework description ****
x_raw = 'HHHTHHTHHHHTTHHHHHHTHHTHHTHHHTHHTTHHHTHHHHHT'  # Raw observation
a = np.array([
    [0.0, 0.5, 0.5], 
    [0.0, 0.8, 0.2],
    [0.0, 0.3, 0.7]
])  # transition matrix a_kl. k=0 is B, k=1 is F, k=2 is L.

e = np.array([
    [0.5, 0.5],
    [0.9, 0.1],
])  # emission probability matrix e_k(b) for states B=0, F=1 and L=2. Column 0 is H, 1 is T 

#****  Pre-processing: replace H with 0 and T with 1 in x ****
x = np.array([ int(s=='T') for s in x_raw ])

#****  Obtain most probable state path: Viterbi algorithm ****
v, pi_m, p_pxpi = viterbi(x, a, e)

#**** Posterior decoding - obtain f and b matrices ****
f, p_f = forward(x, a, e)
b, p_b = backward(x, a, e)
assert np.abs(p_f-p_b)/p_f < 1e-9
p = (p_f + p_b) / 2  # Average probability
p_post = f*b/p  # Posterior probability
# Extract highest probability sequence - not necessarily allowed transitions-wise!
pi_hat = np.argmax(p_post, axis=1)

#**** Compare Viterbi and posterior ****
pp = pprint.PrettyPrinter(indent=4)
print('****** Observed sequence: x ******')
print(f'x: {x_raw}')
print('****** Viterbi: pi_m ******')
pi_m_str = ''.join(translate_to_state(pi_m))
print(r'pi_m: ' + pi_m_str)
print('****** Posterior: pi_hat ******')
pi_hat_str = ''.join(translate_to_state(pi_hat))
print(r'pi_hat: ' + pi_hat_str)
print('****** Fractional overlap Viterbi/Posterior ******')
print(f'Overlap: {(np.sum(pi_m==pi_hat)/len(pi_m)):.4f}')
print('****** Sum of posterior probabilities Sum_k( p(xi=k|x) ) ******')
pp.pprint(np.sum(p_post, axis=1))
# print('****** Posterior probaility p(xi=k|x) ******')
# pp.pprint(p_post)
print('****** First six columns of v ******')
pp.pprint( v[:6,:] )
print('****** P(pi_i=L|x) for the first six observations HHHTHH ******')
pp.pprint( p_post[:6, 2] )
# Plot 
fig, ax = plt.subplots(figsize=(8,6))
ax.plot(pi_m, linewidth=3, linestyle='--', label=r'$\pi^*$, Viterbi')
ax.plot(pi_hat, linewidth=2, linestyle='-', alpha=0.7, label=r'$\hat{\pi}$, Posterior')
# ax.plot(1+p_post[:,1], linewidth=3, color='b', linestyle=':', alpha=0.5, label=r'$P(x_i = F | x$')
# ax.plot(1+p_post[:,2], linewidth=3, color='r', linestyle=':', alpha=0.5, label=r'$P(x_i = L | x$')
ax.grid()
ax.legend(loc='best')
ax.set_xlabel(r'$i$, sequence index')
ax.set_ylabel(r'Predicted HMM state')
plt.locator_params(axis='y', nbins=2)
labels = [item.get_text() for item in ax.get_yticklabels()]
labels[1] = 'F'
labels[2] = 'L'
ax.set_yticklabels(labels)
plt.tight_layout()
plt.savefig('viterbi_post_comparison.png')