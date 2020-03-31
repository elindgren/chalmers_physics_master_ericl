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
A training example based on the "Occasionally dishones casino" example 
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
    e_hf = 0.5          # p head fair coin
    e_hl = 0.1          # p head loaded coin
    e = [e_hl, e_hf]

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
    for i in tqdm(range(1, len(x)+1)):
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
        Given a sequence x, calculate the most probable state sequence pi.
        Perform calculations in log space

        Avoiding underflow: http://wittawat.com/posts/log-sum_exp_underflow.html
    '''
    # Initialization
    f = np.zeros((len(x)+1, 3)) 
    o = 0  # Almost 0, will be used as 0 here to avoid log errors
    #          0  l  f
    f[0, :] = [1, o, o]  # p of most probable path ending in each state, if it ends now
    a_kl = np.array([ [o, a_0k[0], a_0k[1]], [o, 1-a[0], a[0]], [o, a[1], 1-a[1]] ])  # [ [a_00, a_0l, a_0f], [a_l0, a_ll, a_lf], [a_20, a_fl, a_ff] ]
    e_lx = np.array([ [o, o], [e[0], 1-e[0]], [e[1], 1-e[1]] ])  # [ [e_0(h), e_0(t)], [e_l(h), e_l(t)], [e_f(h), e_f(t)] ]
    # Convert to log space
    f[0,:] = np.log(f[0,:])
    a_kl = np.log(a_kl)
    e_lx = np.log(e_lx)
    # Recursion
    for i in tqdm(range(1, len(x)+1)):
        xi = x[i-1]
        for l in range(f.shape[1]):
            # l is the new state, either 0 (loaded) or 1 (fair). 0 is inaccessible p=0
            e_l = e_lx[l, :]
            f_a_kl = f[i-1,:] + a_kl[:, l]
            c = np.max(f_a_kl)
            print(np.sum(np.exp( f_a_kl )))
            f[i, l] =  e_l[xi] + c + np.log(np.sum(np.exp( f_a_kl - c)))
        print(f[i,:])
        print('_____________')
    # Termination
    p_x = np.exp( np.sum( f[-1,:] ) )
    # Traceback
    print(p_x)
    return p_x


# Generate sequence of coin tosses
L = 2
x, pi, e, a, a_0k = casino(L)

# Obtain most probable state sequence - Viterbi algorithm
pi_m = viterbi(x, e, a, a_0k)

# Plot
fig, ax = plt.subplots(figsize=(8,6))
ax.plot(pi_m[1:], c='C1', linewidth=2, alpha=0.7, label=r'$\pi^*$, Viterbi')
# ax.fill_between(np.arange(L), 1, pi_m)
ax.plot(pi, c='k', linestyle='--', linewidth=1, label=r'$\pi$, True sequence')
ax.grid()
ax.legend(loc='best')
ax.set_title('Dishonest casino: Viterbi algorithm')
ax.set_xlabel(r'$i$')
ax.set_ylabel(r'Heads/tails')
plt.savefig(f'viterbi_L={L}.png')

# Obtain probability of sequence - Forward algorithm
p_x = forward(x, e, a, a_0k)

# Print results
print()
print('--------- Results ---------')
print(f'**** Viterbi accuracy: {(np.sum(pi_m == pi) / L):.2f}')
print(f'**** Forward: P(x): {p_x}')