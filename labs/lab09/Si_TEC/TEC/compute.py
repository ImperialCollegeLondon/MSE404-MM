import numpy as np
import matplotlib.pyplot as plt

from tec import *

T = np.arange(50, 2000, 300)
enfile = "etot_v_vol.dat"
freqfiles = np.array(["Si00-q20.freq",\
             "Si01-q20.freq",\
             "Si02-q20.freq",\
             "Si03-q20.freq",\
             "Si04-q20.freq",\
             "Si05-q20.freq"])

V = alpha(enfile, freqfiles, T)
## Make a plot
plt.plot(T, V, color="k", label="Vvs.T")
plt.xlabel(r"Temperature (K)", fontsize=16)
plt.ylabel(r"Volume(Ang^3)", fontsize=16)
plt.savefig("VvsT.png", dpi=300)
plt.show()
