# External imports
import numpy as np
from pandas import *

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


# Sequences
x_s = 'ACGGGTAATTAAGCCGA'
y_s = 'CCGATATTTTAGCA'
s = np.array([
    [2, -2, 1, 1],
    [-2, 12, -3, -2],
    [1, -3, 5, 0],
    [1, -2, 0, 3]
])
x = translate_alphabet(x_s)
y = translate_alphabet(y_s)

# Setup
d = 2
F = np.zeros((len(x), len(y)))
T = np.zeros((len(x), len(y)))  # traceback table

# Initialization
for i in range(len(x)):
    F[i,0] = -d*i
for j in range(len(y)):
    F[0,j] = -d*j

# Algorithm
for i in range(1, len(x)):
    for j in range(1, len(y)):
        cases = [ F[i-1, j-1]+s[x[i], y[j]], F[i-1, j]-d, F[i, j-1]-d ] 
        choice = np.argmax(cases) # 0 = i-1,j-1, 1=i-1,j (xi aligned to gap), 2=i,j-1 (yi aligned to gap)
        F[i,j] = cases[choice] # i is the column index, j is the row index
        T[i,j] = choice

# Traceback
x_al = []
y_al = []
path = []
i = len(x)-1
j = len(y)-1
while i>-1 or j>-1:
    path.append([i,j])
    if i==len(x)-1 and j==len(y)-1:
        c=0 # The last characters must be aligned
    else:
        c=T[i,j]
    if c==0:
        # xi and yi aligned
        x_al.append(x_s[i])
        y_al.append(y_s[j])
        i -= 1
        j -= 1
    elif c==1:
        # xi aligned to gap in y
        x_al.append(x_s[i])
        y_al.append('-')
        i -= 1
        j = j # no change in y, gap
    elif c==2:
        # yi aligned to gap in x
        x_al.append('-')
        y_al.append(y_s[j])
        i = i # no change in x - gap in 
        j -= 1 
x_al = ''.join(list(reversed(x_al)))
y_al = ''.join(list(reversed(y_al)))

# Convert F to string and mark path
Fs = []
for i in range(len(x)):
    row = []
    for j in range(len(y)):
        if [i,j] in path:
            row.append(f'[{int(F[i,j])}]')
        else:
            row.append(f' {int(F[i,j])} ')
    Fs.append(row)
Fd = DataFrame(Fs)
Fd.index = [xi for xi in x_s]
Fd.columns = [yi for yi in y_s]
print(Fd)


print(f'Optimal global alignment with score {F[-1,-1]}:')
print(f'x: {x_al}')
print(f'y: {y_al}')