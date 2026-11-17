import sys
import numpy as np

filename = sys.argv[1] 
T = float(sys.argv[2])

data = np.loadtxt(filename, skiprows=1)
freqs = data[:,0]
dos = data[:,1]
constant1 = ???
constant2 = ???
dw = freqs[1] - freqs[0] 

total_energy = 0.0
for w, rho in zip(freqs, dos):

  if abs(w) < 1e-5:
    continue 

  bose = 1/(np.exp(constant1*w/T) - 1)
  total_energy += constant2*w*rho*(0.5 + bose)*dw

print(total_energy)
