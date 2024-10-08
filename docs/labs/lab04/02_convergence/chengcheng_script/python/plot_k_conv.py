import matplotlib.pyplot as plt
import numpy as np

# read-in data
data = np.loadtxt('etot_v_nkpt.dat')
# the data is structured as
# num_kpt | etot [Ry]

# plot the convergence  
plt.plot(data[:,0], data[:,1], marker='o')

# set xlabel
plt.xlabel('number of k-points along each direction')
plt.ylabel('Total energy [Ry]')

plt.savefig('k_conv.pdf')
