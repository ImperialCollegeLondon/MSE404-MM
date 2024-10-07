#Materials Modelling course
#--------------------------

import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.interpolate import CubicSpline
from thermo import *


def Ry2mev():
  return 13.60569301*1000


def energy(enfile):
  """
  From a certain file extracts the total energy
  as a function of volume;
  """
  # Energy is in Rydberg
  V, E = np.loadtxt(enfile, unpack=True)
  return V, E*Ry2mev()


def phonon_contribution(freqfiles, T):
  """
  For the set of phonon frequecy files compte the
  phonon energies and entropy (in units of meV); 
  Be mindful without the unit-conversion
  """
  E_ph_V = np.zeros((T.shape[0], freqfiles.shape[0]))
  TS_V = np.zeros((T.shape[0], freqfiles.shape[0]))
  for i in range(freqfiles.shape[0]):
    print("Phonon file:%s"%(freqfiles[i]))
    for j in range(T.shape[0]):
      Thermo = Thermodynamic(freqfiles[i], T[j])
      E_ph_V[j][i] = Thermo.get_E_T()
      TS_V[j][i] = Thermo.get_TS()
  return E_ph_V, TS_V 


def minimize_H(E, E_ph, TS, V):
  """
  Computes Helmohltz free energy at a specific volume
  for a given temperature.
  @input
    E: Energy for a series of volume; Temperature independent
    E_ph: Phonon energy at that temperature for a series of volumes
    S: Entropy at that same temperature
    V: Volumes
  """
  H = E + E_ph - TS
  # Interpolate to get the minimum volume
  V_ = np.arange(np.min(V), np.max(V), 0.001)
  H_interp = CubicSpline(V, H)
  H_ = H_interp(V_)
  V_min = V_[np.where(H_ == H_.min())]
  return V_min


def thermal_expansion_coeff(V, E_V, E_ph_V, TS_V, T):
  """
  Computes thermal expansion coefficient
  @input:
    E_V: Energy for each volume
    E_ph_V: Phonon energy for a series of temperature
            at each volume.
            Dimensionally, N_TxN_V
    S_V: Entropy for a series of temperature at each volume
  """
  V_opt = np.zeros((E_ph_V.shape[0]))
  for i in range(E_ph_V.shape[0]):
    V_opt[i] = minimize_H(E_V, E_ph_V[i], TS_V[i], V)
  return V_opt 


def alpha(enfile, freqfiles, T):
  """
  Computes thermal expansion coefficients at T
  """
  V, E_V = energy(enfile)
  E_ph_V, TS_V = phonon_contribution(freqfiles, T)
  V_opt = thermal_expansion_coeff(V, E_V, E_ph_V, TS_V, T)
  return V_opt
