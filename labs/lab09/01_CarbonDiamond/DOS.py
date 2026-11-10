import numpy as np
import matplotlib.pyplot as plt
dos = np.loadtxt("CD.dos", skiprows=1)
plt.plot(dos[:,0], dos[:,1])
plt.xlabel(r'Frequency ($\mathrm{cm}^{-1}$)',fontsize=12)
plt.ylabel(r'DoS (states/$\mathrm{cm}^{-1}$)',fontsize=12)
plt.savefig("CD_dos.png")
plt.show()
