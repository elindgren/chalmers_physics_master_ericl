# Local imports

# Global imports
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

# Set plot params
plt.rc('font', size=18)          # controls default text sizes
plt.rc('axes', titlesize=18)     # fontsize of the axes title
plt.rc('axes', labelsize=18)     # fontsize of the x and y labels
plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
plt.rc('ytick', labelsize=18)    # fontsize of the tick labels
plt.rc('legend', fontsize=18)    # legend fontsize

'''
A training example based on the "Occasionally dishonest casino" example 
from  chapter 3 of Biological Sequence analysis, focusing on HMMs.
The specific example is random tosses of coins, and the task is to find 
    - most likely state sequence
    - probability of a certain sequence
    - most probable state for a certain observation
    - estimates for HMM model parameters, in two cases
        * Known state sequence => MLE estimates
        * Unknown state sequence: EM with Baum-Welch
'''


def casino(N):
    '''
        Generates a sequence lenght N of random contosses, with either a weighted coin or with a fair one. 
        Define: Head = h = 0, Tails = t = 0
    '''
    np.random.seed(1)
    
    # Data structures
    coin_seq = np.zeros(N, dtype=int)
    state_seq = np.zeros(N, dtype=int)

    # Transition probabilities 
    a_lf = 0.1         # p loaded => fair
    a_fl = 0.01         # p fair => loaded
    a = [a_lf, a_fl]
    
    a_0l = 0.3  # p start
    a_0f = 1 - a_0l
    a_0k = [a_0l, a_0f]
    state = np.random.choice(a=[0,1], size=1, p=a_0k)[0]  # 0 = loaded, 1 = fair, 

    # Emission probabilities
    e_fh = 0.5          # p head fair coin
    e_lh = 0.1          # p head loaded coin
    e = [e_lh, e_fh]

    for i in range(N):
        # go to new state based on previous state
        transition = np.random.choice([True, False], size=1, p=[a[state], 1-a[state]])[0]
        if transition:
            state = not state
        # emit value 
        emission = np.random.choice([0,1], size=1, p=[e[state], 1-e[state]])[0]
        # save emitted value and state
        coin_seq[i] = emission
        state_seq[i] = state
    return coin_seq, state_seq, e, a, a_0k


def viterbi(x, e, a, a_0k):
    ''' 
        Given a sequence x, calculate the most probable state sequence pi.
        Perform calculations in log space
    '''
    # Initialization
    pi = np.zeros(len(x))
    v = np.zeros((len(x)+1, 3)) 
    #          0  l  f
    v[0, :] = [1, 0, 0]  # p of most probable path ending in each state, if it ends now
    a_kl = np.array([ [0, a_0k[0], a_0k[1]], [0, 1-a[0], a[0]], [0, a[1], 1-a[1]] ])  # [ [a_00, a_0l, a_0f], [a_l0, a_ll, a_lf], [a_20, a_fl, a_ff] ]
    e_lx = np.array([ [0, 0], [e[0], 1-e[0]], [e[1], 1-e[1]] ])  # [ [e_0(h), e_0(t)], [e_l(h), e_l(t)], [e_f(h), e_f(t)] ]
    # Convert to log space
    v[0,:] = np.log(v[0,:])
    a_kl = np.log(a_kl)
    e_lx = np.log(e_lx)
    # Recursion
    for i in range(1, len(x)+1):
        xi = x[i-1]
        for l in range(v.shape[1]):
            # l is the new state, either 0 (loaded) or 1 (fair). 0 is inaccessible p=0
            e_l = e_lx[l, :]
            v_a_kl = v[i-1,:] + a_kl[:, l]
            v[i, l] =  e_l[xi] + np.max( v_a_kl ) # sufficient to max since log is monotonous
            pi[i-1] = np.argmax( v_a_kl ) - 1 # remove the extra state
    # Termination
    p_xpi = np.exp(np.max(v[-1,:]))
    # Traceback
    return pi

def forward(x, e, a, a_0k):
    ''' 
        Given a sequence x, calculate its probability p(x).
        Perform calculations in log space

        Avoiding underflow: http://wittawat.com/posts/log-sum_exp_underflow.html
    '''
    # Initialization
    f = np.zeros((len(x)+1, 3))  # p of most probable path ending in each state, if it ends now
    #          0  l  f
    f[0, :] = [1, 0, 0] 
    a_kl = np.array([ [0, a_0k[0], a_0k[1]], [0, 1-a[0], a[0]], [0, a[1], 1-a[1]] ])  # [ [a_00, a_0l, a_0f], [a_l0, a_ll, a_lf], [a_f0, a_fl, a_ff] ]
    e_lx = np.array([ [0, 0], [e[0], 1-e[0]], [e[1], 1-e[1]] ])  # [ [e_0(h), e_0(t)], [e_l(h), e_l(t)], [e_f(h), e_f(t)] ]
    # Convert to log space
    f[0,:] = np.log(f[0,:])
    a_kl = np.log(a_kl)
    e_lx = np.log(e_lx)
    # Recursion
    for i in range(1, len(x)+1):
        xi = x[i-1]
        for l in range(f.shape[1]):
            # l is the new state, either 0 (loaded) or 1 (fair). 0 is inaccessible p=0
            e_l = e_lx[l, :]
            f_a_kl = f[i-1,:] + a_kl[:, l]
            c = np.max(f_a_kl)
            fa_sum = c + np.log(np.sum(np.exp( f_a_kl - c)))
            # print( a_kl[l, :])
            if e_l[xi] == -np.inf or fa_sum == -np.inf :
                f[i, l] = -np.inf
            else:
                f[i, l] =  e_l[xi] + fa_sum
    # Termination
    a_k0 = np.array([0, 0.5, 0.5])
    p_x = np.sum( np.exp(f[-1,:]) * a_k0 )
    # p_x = np.sum(np.exp( f[-1, :] ))
    print(p_x)
    # Unlog
    f = np.exp(f)
    return p_x, f


def backwards(x, e, a, a_0k):
    ''' 
        Given a sequence x, calculate b in order to calculate the probability that 
        the HMM is in state k at char i p( pi_i = k | x ).
        Perform calculations in log space
    '''
    # Initialization
    b = np.zeros((len(x)+1, 3)) 
    #          0  l  f
    b[-1, :] = [0, 0.5, 0.5]  # bk(L) = a_k0, i.e.prob to go to ending state
    a_kl = np.array([ [0, a_0k[0], a_0k[1]], [0, 1-a[0], a[0]], [0, a[1], 1-a[1]] ])  # [ [a_00, a_0l, a_0f], [a_l0, a_ll, a_lf], [a_f0, a_fl, a_ff] ]
    e_lx = np.array([ [0, 0], [e[0], 1-e[0]], [e[1], 1-e[1]] ])  # [ [e_0(h), e_0(t)], [e_l(h), e_l(t)], [e_f(h), e_f(t)] ]
    # Convert to log space
    b[0,:] = np.log(b[0,:])
    a_kl = np.log(a_kl)
    e_lx = np.log(e_lx)
    # Recursion
    for i in reversed(range(0, b.shape[0]-1)):
        xi = x[i]
        for l in range(b.shape[1]):
            # l is the new state, either 0 (loaded) or 1 (fair). 0 is inaccessible p=0
            e_l = e_lx[l, :]
            b_a_kl = b[i+1,:] + a_kl[:, l]
            c = np.max(b_a_kl)
            ba_sum = c + np.log(np.sum(np.exp( b_a_kl - c)))
            # print( a_kl[l, :])
            if e_l[xi] == -np.inf or ba_sum == -np.inf:
                b[i, l] = -np.inf
            else:
                b[i, l] =  e_l[xi] + ba_sum
    # Termination
    x1 = x[0]  # stupid algorithm count from 1 >.<
    a_0l = np.array( [0, a_0k[0], a_0k[1]] )
    p_x = 0
    for l in range(b.shape[1]):
        e_l = e_lx[l, :]
        p_x +=  a_0l[l] * np.exp( e_l[x1] + b[0,l] ) 
    print(p_x)
    # Unlog
    b = np.exp(b)
    return p_x, b


def baum_welch(X):
    '''
        Compute estimates for all HMM model parameters using EM algorithm.
        X is a matrix containing all {x^j}
    '''
    # Intialisation - random model parameters
    a = np.array( [np.random.random_sample(), np.random.random_sample()] )
    e = np.array( [np.random.random_sample(), np.random.random_sample()] )
    a_0l = np.random.random_sample()
    a_0k = [0, a_0l, 1-a_0l]
    # Recurrence - Iterate until convergence
    for s in range(1):
        print(f'---------- Iteration {s} ----------')
        A = np.zeros((3,3))  # Expectation values of the various transitions - size k x l
        E = np.zeros((3,2))

        a_kl = np.array([ [0, a_0k[0], a_0k[1]], [0, 1-a[0], a[0]], [0, a[1], 1-a[1]] ])  # [ [a_00, a_0l, a_0f], [a_l0, a_ll, a_lf], [a_f0, a_fl, a_ff] ]
        e_lx = np.array([ [0, 0], [e[0], 1-e[0]], [e[1], 1-e[1]] ])
        for j in range(X.shape[0]):
            # Get f and b
            x_j = X[j,:]
            p_xf, f = forward(x, e, a, a_0k)
            p_xb, b = backwards(x, e, a, a_0k)
            
            # Log parameters
            a_kl = np.log(a_kl)
            e_lx = np.log(e_lx)
            f = np.log(f)
            b = np.log(b)
            p_xf = np.log(p_xf)
            p_xb = np.log(p_xb)
            for k in range(A.shape[0]):
                # Estimate A
                for l in range(A.shape[0]):
                    e_l = e_lx[l, :]
                    for i in range(len(x)-1):
                        e_lxi1 = e_l[x_j[i+1]]
                        fb_sum = f[i,k] + b[i+1, l] - p_xf
                        print(f[i,k])
                        if e_lxi1 == -np.inf or fb_sum == -np.inf :
                            A[k, l] = np.exp( -np.inf )
                        else:
                            A[k,l] += np.exp( e_lxi1 + fb_sum )
                # Estimate E
                for s in range(1):
                    s = int(s)
                    # Find all s=h or s=f
                    x_idx = np.where(x == s)[0]
                    for i in x_idx:
                        E[k, s] += np.exp( f[i,k] + b[i,k] - p_xf)
        # Calculate new parameter estimates - MLE estimation
        a_kl = np.zeros((3,3))
        e_kb = np.zeros((3,2))
        for k in range(3):
            for l in range(3):
                a_kl[k,l] = A[k,l] / np.sum( A[k,:] )
            for s in range(1):
                s = int(s)
                e_kb[k,s] = E[k,s] / np.sum( E[k,:] )
            # Normalize to probabilities
            a_kl[k,:] /= np.sum( a_kl[k,:] ) # sum_l a_kl = 1, must go somewhere
            e_kb[k,:] /= np.sum( e_kb[k,:] )
        
        print(a_kl)
        
        a = np.array( [a_kl[1,2], a_kl[2,1]] )
        e = np.array( [e_kb[1,0], e_kb[2,0]] )
        a_0k = [0, a_kl[0,1], a_kl[0,2]]
        # Calculate Likelihood of model
    return None


# Generate sequence of coin tosses
L = 300
x, pi, e, a, a_0k = casino(L)

# Obtain most probable state sequence - Viterbi algorithm
pi_m = viterbi(x, e, a, a_0k)


# Obtain probability of sequence - Forward algorithm
p_xf, f = forward(x, e, a, a_0k)

# Obtain probability of state for xi - Backwards algorithm
p_xb, b = backwards(x, e, a, a_0k)
p_pii = f * b / p_xb
# Error somewhere, rescale p_pii so sums to 1
# p_pii /= np.sum(p_pii[1,:])

# HMM Parameter Estimation - Baum-Welch
s = baum_welch(np.array([x]))


# Print results
print()
print('--------- Results ---------')
print(f'**** Viterbi accuracy: {(np.sum(pi_m == pi) / L):.2f}')
print(f'**** Forward: P(x): {p_xf}')
print(f'**** Backward: P(pi_100 = k | x): {p_pii[1,:]}')

# Plot
fig, ax = plt.subplots(figsize=(8,6))

ax.plot(pi, c='k', linestyle='--', linewidth=1, label=r'$\pi$, True sequence')
ax.plot(pi_m[1:], c='C1', linewidth=2, alpha=0.7, label=r'$\pi^*$, Viterbi')

ax.plot(p_pii[:,1], c='r', linestyle=':', linewidth=3, alpha=0.7, label=r'$p(x, \pi_i=l)$, Posterior')
ax.plot(p_pii[:,2], c='c', linestyle=':', linewidth=2, alpha=0.7, label=r'$p(x, \pi_i=f)$, Posterior')

ax.grid()
ax.legend(loc='lower left')
ax.set_title('The Dishonest casino')
ax.set_xlabel(r'$i$')
ax.set_ylabel(r'Heads/tails')
plt.tight_layout()
plt.savefig(f'viterbi_L={L}.png')
# plt.show()