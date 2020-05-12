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
x_s = 'CGATTA'            #! change
y_s = 'ATTTC'               #! change
s = np.array([
    [1, -2, -2, -2],
    [-2, 2, -2, -2],
    [-2, -2, 2, -2],
    [-2, -2, -2, 1]
])                          #! change
x = translate_alphabet(x_s)
y = translate_alphabet(y_s)

# Setup
d = 4
F = np.zeros((len(x)+1, len(y)+1))
T = np.zeros((len(x)+1, len(y)+1))  # traceback table

# Initialization
for i in range(F.shape[0]):
    F[i,0] = 0
for j in range(F.shape[1]):
    F[0,j] = 0

# Algorithm
for i in range(1, F.shape[0]):
    for j in range(1, F.shape[1]):
        cases = [ 0, F[i-1, j-1]+s[x[i-1], y[j-1]], F[i-1, j]-d, F[i, j-1]-d ] 
        choice = np.argmax(cases) # 0 = i-1,j-1, 1=i-1,j (xi aligned to gap), 2=i,j-1 (yi aligned to gap)
        F[i,j] = cases[choice] # i is the column index, j is the row index
        T[i,j] = choice

# Traceback
x_al = []
y_al = []
path = []
max_idx = np.where(F == F.max())
i = max_idx[0][0]
j = max_idx[1][0]
while i>0 or j>0:
    c=T[i,j]
    path.append([i,j])
    if c==0:
        # the case of Fij = 0, i.e. beginning of alignment. End
        # x_al.append(x_s[i-1])
        # y_al.append(y_s[j-1])
        break
    elif c==1:
        # xi and yi aligned
        x_al.append(x_s[i-1])
        y_al.append(y_s[j-1])
        i -= 1
        j -= 1
    elif c==2:
        # xi aligned to gap in y
        x_al.append(x_s[i-1])
        y_al.append('-')
        i -= 1
        j = j # no change in y, gap
    elif c==3:
        # yi aligned to gap in x
        x_al.append('-')
        y_al.append(y_s[j-1])
        i = i # no change in x - gap in 
        j -= 1
x_al = ''.join(list(reversed(x_al)))
y_al = ''.join(list(reversed(y_al)))

# Convert F to string and mark path
Fs = []
for i in range(F.shape[0]):
    row = []
    for j in range(F.shape[1]):
        if [i,j] in path:
            row.append(f'[{int(F[i,j])}]')
        else:
            row.append(f' {int(F[i,j])} ')
    Fs.append(row)
Fd = DataFrame(Fs)
rows = ['-']
tmp = [rows.append(xi) for xi in x_s]
Fd.index = rows
cols = ['-']
tmp = [cols.append(yi) for yi in y_s]
Fd.columns = cols
print(Fd)


print(f'Optimal local alignment with score {F.max()}:')
print(f'x: {x_al}')
print(f'y: {y_al}')