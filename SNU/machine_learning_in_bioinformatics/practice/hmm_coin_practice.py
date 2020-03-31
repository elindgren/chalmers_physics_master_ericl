# Local imports

# Global imports
import numpy as np
import matplotlib.pyplot as plt

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
    coin_seq = np.zeros(len(N))
    state_seq = np.zeros(len(N))

    # Transition probabilities 
    a_fl = 0.90         # p loaded => fair
    a_lf = 0.05         # p fair => loaded
    a = [a_fl, a_lf]
    
    # Emission probabilities
    e_hf = 0.5          # p head fair coin
    e_hl = 0.1          # p head loaded coin
    e = [e_hl, e_hf]

    state = np.random.chose(a=[0,1], size=1, p=[0.5, 0.5])  # chose either state, with 50/50 chance. 1 = fair, 0 = loaded
    for i in range(N):
        # go to new state based on previous state
        transition = np.random.chose([True, False], size=1, p=[a[state], 1-a[state]])
        if transition:
            state = not state
        # emit value 
        emission = np.random.chose([0,1], size=1, p=[e[state], 1-e[state]])
        # save emitted value and state
        coin_seq[i] = emission
        state_seq[i] = state


    return coin_seq, state_seq


# Generate sequence of coin tosses
coin_seq, state_seq = casino(10)

print(coin_seq)
print(state_seq)