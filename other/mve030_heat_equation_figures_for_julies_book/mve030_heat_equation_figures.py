#  |\__/,|   (`\
#  |_ _  |.--.) )
#  ( T   )     /
# (((^_(((/(((_/
#
# Thank you for checking out this script!
# Author: Eric Lindgren, F16, 2021

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# Set plot params
plt.rc('font', size=16)          # controls default text sizes
plt.rc('axes', titlesize=16)     # fontsize of the axes title
plt.rc('axes', labelsize=18)     # fontsize of the x and y labels
plt.rc('xtick', labelsize=16)    # fontsize of the tick labels
plt.rc('ytick', labelsize=16)    # fontsize of the tick labels
plt.rc('legend', fontsize=16)    # legend fontsize

def alpha0(x,f):
    return np.trapz(y=f, x=x) / (2*np.pi)

def alpha(n,x,f):
    return np.trapz(y=f*np.cos(n*x), x=x) / np.trapz(y=np.abs(np.cos(n*x))**2, x=x)

def beta(n,x,f):
    return np.trapz(y=f*np.sin(n*x), x=x) / np.trapz(y=np.abs(np.sin(n*x))**2, x=x)

def un(x,t,n,k,f):
    if n == 0:
        alpha_n = alpha0(x,f)
        beta_n = 0
    else: 
        alpha_n = alpha(n,x,f)
        beta_n = beta(n,x,f)
    return np.exp(-n**2 * t*k) * (alpha_n * np.cos(n*x) + beta_n * np.sin(n*x))

def u(x,t,n_max,k,f):
    u_xt = np.zeros(len(x))
    for n in range(n_max):
        u_xt += un(x,t,n,k,f)
    return u_xt


def circular_rod_heat_solution(n_max, x, f, flabel, filename):
    '''
    Plots the solution to the heat equation on a circular one dimensional rod, for a given vector of x-values [-pi, pi]
    and an initial function f defined on those x-values with a corresponding label flabel. 
    Args:
      n_max (number): Integer for the maximum n (i.e. the number of basis functions) that will be used.
      x (list): x-domain, vector of values in the range [-pi, pi].
      f (list): an initial value function defined on the x-domain.
      flabel (string): Label of the initial value function for naming the output file.
    '''
    # Define constants - here one can control the number of ns, i.e. the accuracy of the solution
    k = 1
    ts = [0,1,100]

    # Calculate heat equation solution for the timesteps
    sols = np.zeros((len(ts), len(x)))

    for i,t in enumerate(ts):
        sol = u(x,t,n_max,k,f)
        sols[i,:] = sol

    # Plot solutions
    fig, ax = plt.subplots(figsize=(8,6))

    cs = ['r', 'y', 'b']
    styles = ['-', '-.', ':']
    ax.plot(x, f, c='k', linestyle='--', linewidth=2, alpha=1, label=r'$f(x)=$'+flabel)
    for i,sol in enumerate(sols):
        ax.plot(x, sol, c=cs[i], alpha=0.6, linestyle=styles[i], linewidth=2, label=r'$u(x,t=$'+f'{ts[i]:.0f}'+r'$)$')
    ax.legend(loc='best')
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$u(x,t)$, arb. units')
    tick_places = [-np.pi, -np.pi/2, 0, np.pi/2, np.pi]
    my_xticks = [r'$-\pi$',r'$-\frac{\pi}{2}$','0',r'$\frac{\pi}{2}$',r'$\pi$']
    plt.xticks(tick_places, my_xticks)

    # Disable y ticks
    # frame1 = plt.gca()
    # frame1.axes.get_yaxis().set_ticks([])

    # Save figure
    plt.tight_layout()
    plt.savefig(f'./images/pdf/f_{filename}_n={n_max}_k={k}.pdf') # Save as pdf
    plt.savefig(f'./images/png/f_{filename}_n={n_max}_k={k}.png') # and png
    plt.close()


def main():
    n_max = 10 # number of ns - i.e. NOT including n = n_max
    x = np.linspace(-np.pi, np.pi, 500)
    # Initial function values - add/modify these and the labels/filenames below appropriately
    fs = [
        2*np.ones(len(x)),         # Constant
        x,                         # Linear
        np.sin(x),                 # Sine
        np.heaviside(x, 0.5),      # Heaviside - x2 is the value at f(x=0)
        x**2,                      # Quadratic
        np.cosh(x),                # Cosh
        np.abs(x),                 # Absolute value
    ]
    flabels = [
        'constant',
        r'$x$',
        r'$\sin(x)$',
        r'$H(x)$',
        r'$x^2$',
        r'$\cosh(x)$',
        r'$|x|$'
    ]
    filenames = [
        'constant',
        'linear',
        'sine',
        'heaviside',
        'quadratic',
        'cosh',
        'abs'
    ]

    for i, f in enumerate(fs):
        circular_rod_heat_solution(n_max, x, f, flabels[i], filenames[i])


if __name__ == "__main__":
    main()