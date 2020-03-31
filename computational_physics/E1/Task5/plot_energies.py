# plot the energies
# Created by Martin Gren 2014-10-25.

# imports
import matplotlib.pylab as plt
import numpy as np

# input file
filename = 'energy.dat'

# import data
data = np.loadtxt(filename)

# size of figure
plt.figure(figsize=(8,6))

# plot
plt.plot(data[:,0], data[:,1],'-',label='Total energy')
plt.plot(data[:,0], data[:,2],'-',label='Potential energy')
plt.plot(data[:,0], data[:,3],'-',label='Kinetic energy')

# labels
plt.xlabel('Time / [dim. unit]', fontsize=12)
plt.ylabel('Energy / [dim. unit]', fontsize=12)

# legend
plt.legend()
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
plt.setp(ltext, fontsize=12) 

# axis limits
plt.xlim([0,50])
plt.ylim([min(data[:,3])-0.002,max(data[:,1])+0.002])

# tick fontsize
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

# save and display the plot
plt.savefig('energies.pdf')
plt.show()
