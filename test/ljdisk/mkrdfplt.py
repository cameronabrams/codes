import matplotlib.pyplot as plt
import numpy as np

dat=np.loadtxt('rdf.dat')
r=dat[:,0]
gr=dat[:,1]

plt.xlim(0,10)
plt.ylim(0,3)
plt.xlabel('r')
plt.ylabel('g(r)')

plt.plot(r,gr,'r-')
plt.savefig('rdf.png')

