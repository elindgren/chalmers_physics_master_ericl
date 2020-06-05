# External imports
import numpy as np

'''
    Calculates the numerical vacancy concentration as per SIA and Vacancy. Do not consider
    Frenkel defects (vacancy + interstitial) - this only works for SIA and Schottky defects.
    Problem: 
    "We have another system composed of 10^23 atoms. The sort of atoms is
    the same with that in Quiz-1 (and thus the vacancy formation energy is 2.0 eV).
    Determine the equilibrium vacancy concentration in this system at 300 K and
    1000 K. You can assume that all the vacancies are formed as Schottky defects."
'''


def c(N, E_def, T):
    kB = 8.6173e-5 # eV/K
    return N*np.exp(-E_def/(kB*T))


N = 1e23 # number of atoms
# Schottky (Vacancy)
E_vac = 1.7 # eV
T1 = 300
T2 = 1000
c_vac_300 = c(N, E_vac, T1)
c_vac_1000 = c(N, E_vac, T2)

print(f'Schottky (vacancy) concentration at {T1} K: {c_vac_300:e} --- {T2} K: {c_vac_1000:e}')

# SIA, 111
E_sia = 4.0 # eV
T1 = 300
T2 = 1000
c_vac_300 = c(N, E_sia, T1)
c_vac_1000 = c(N, E_sia, T2)

print(f'SIA111 concentration at {T1} K: {c_vac_300:e} --- {T2} K: {c_vac_1000:e}')