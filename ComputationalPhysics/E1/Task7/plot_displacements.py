# plot the displacements
# Created by Martin Gren 2014-10-25.

# imports
import matplotlib.pylab as plt
import numpy as np

# input file
filename = 'disp.dat'

# import data
data = np.loadtxt(filename)

# size of figure
plt.figure(figsize=(8,6))

# plot
plt.plot(data[:,0], data[:,1],'-',label='Atom 1')
plt.plot(data[:,0], data[:,2],'-',label='Atom 2')
plt.plot(data[:,0], data[:,3],'-',label='Atom 3')

# labels
plt.xlabel('Time / [dim. unit]', fontsize=12)
plt.ylabel('Displacement / [dim. unit]', fontsize=12)

# legend
plt.legend()
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
plt.setp(ltext, fontsize=12) 

# axis limits
plt.xlim([0,50])

# tick fontsize
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

# save and display the plot
plt.savefig('displacements.pdf')
plt.show()
