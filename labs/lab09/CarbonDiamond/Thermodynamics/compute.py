# Materials Modelling course
#--------------------------

from thermo import *
import numpy as np
import matplotlib.pyplot as plt

# Thermo class
# 1st input: Your band-structure file for a fine-grid
# 2nd input: temperature in K
Thermo = Thermodynamic("CD-fine.freq", 20)

# Thermal properties at specific T
#================================

# Total Energy calculations
#-------------------------
E_T = Thermo.get_E_T()
print("Total energy: %.6f meV"%(E_T))

# Specific heat calculations
#--------------------------
Cv_T = Thermo.get_Cv_T()
print("Constnat volume specific heat: %.6f"%(Cv_T))

# Helmohtz-free energy calculations
#---------------------------------
H_T = Thermo.get_H_T()
print("Helmohtz free energy: %.6f meV"%(H_T))

# Entropy calculations
#---------------------
S_T = Thermo.get_S_T()
print("Entropy=%.6f"%(S_T))


# Temperature dependent thermal properties
#========================================
# Implement Temperature dependence 

