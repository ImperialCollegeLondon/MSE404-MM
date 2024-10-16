# Materials modelling course
# Author: Indrajit Maity
# Email: i.maity@imperial.ac.uk
#------------------------------


# Import all the things we need
import numpy as np
import sys


# This is a clss
class Thermodynamic(object):
  """
  Thermodynamic class
  Units: Always in meV 
  """

  def __init__(self, inp_f, T):
    """
  @input
    inp_f: Input file containing the bands;
           This is the same output file of your 
           density-of-state calculations
    T: Temperature, at which you want to calculate things
    """
    self.inp_f = inp_f
    # Kelvin to meV conversion
    self.T = T*self.temp2mev()
    print("Temperature :%d K"%(self.T/self.temp2mev()))
    # set omega_nq and q-grid for rest 
    self.omega_nq, self.qgrid = self.get_omega_nq()
        


  def temp2mev(self):
    """
    Scale factor to convert temperature to meV
    With this factor, 300 K = 25.85 meV
    """
    return 0.08617328149740566


  def set_nbnd(self, string):
    """
    Extracts integer number, nbnd from a given line
    @input
      string from which nbnd to be determined
    @output
      nbnd: Integer
    """
    # First separate the string with commas
    s = string.split(",")
    # Find "nbnd=" string
    for i in range(len(s)):
      if "nbnd=" in s[i]:
        nbnd= int(s[i].split("=")[1])
    return nbnd

  
  def set_nks(self, string):
    """
    Extracts integer number, nks from a given line
    @input
      string from which nks to be determined
    @output
      nks: Integer
    """
    # First separate the string with commas
    s = string.split(",")
    # Find "nbnd=" string
    for i in range(len(s)):
      if "nks=" in s[i]:
        nks = int(s[i].split("=")[1])
    return nks


  def get_omega_nq(self):
    """
    Extracts the phonon bands, \omega_nq , computed
    at a dense q-grid. 
    @output
      ->Phonon frequencies array, omega_nq
      Frequencies are in meV
      ->qgrid, momentum grid.
    """
    f = open(self.inp_f, "r")
    lines = f.readlines() 
    f.close()

    # Extract number of bands and nkpoints    
    nbnd = self.set_nbnd(lines[0].replace("/",""))
    nkp = self.set_nks(lines[0].replace("/",""))

    # Number of lines for each k-points
    # 10 is standard and thus, hard-coded
    num_l = 10
    if nbnd > num_l:
      l = nbnd%num_l
      if l != 0.0:
        nline = nbnd//num_l + 1
      else:
        nline = nbnd//num_l
    else:
      nline = 1

    # Initialize phonon bands and momentum vectors
    omega_nq = np.zeros((nkp, nbnd))
    qgrid = np.zeros((nkp, 3))

    # For all the lines except 1st-line
    for i in range(1, len(lines), nline+1):
      # Loop for k-points and bands
      for k in range(i, i+nline+1):
        # 1st-point is always momentum, belongs to qgrid
        if k == i:
          ind = k//(nline+1)
          for l in range(3):
            qgrid[ind][l] = eval(lines[k].split()[l])
        # rest of the points are bands
        else:
          for p in range(len(lines[k].split())):
            omega_nq[ind][p]=eval(lines[k].split()[p])
    # Return the omega_nq and qgrid
    return np.moveaxis(omega_nq, 0, -1)*self.cminv2mev(), qgrid
 

  def befactor(self, omega):
    """
    @output
      Bose-Einstein distribution function for 
      a particular phonon mode at a particular
      T;
    """
    if omega < self.T:
      # For very small phonon frequencies, such as 
      # acoustic phonon modes, set omega = 0.01 meV.
      if omega <= 10**-4:
        omega = 0.001
        #print("NOTE: Very-low frequency phonon modes are replaced\n by a non-zero frequency")
    return 1/(np.exp(omega/self.T) - 1)


  def cminv2mev(self):
    """
    Conversion factor from cm-1 to meV;
    1 meV = 8.065610 cm-1 
    """
    return 0.123983


  def get_E_T(self):
    """
    Computes total energy at a Temperature
    @output
      E_T, in units of meV
      NOTE: E_T/nkp is returned
    """
    E_T = 0.0
    # For every band n
    for n in range(self.omega_nq.shape[0]):
      # For every band q
      for q in range(self.omega_nq.shape[1]):
        omega = self.omega_nq[n][q]
        E_T = E_T + \
              (omega*\
              (0.5 + self.befactor(omega)))
    return E_T/self.omega_nq.shape[1]
 

  def get_H_T(self):
    """
    Get Helmohtz free energy at T
    @output
      H_T
      NOTE: H_T/nkp is returned
    """
    H_T = 0.0
    # For every band n
    for n in range(self.omega_nq.shape[0]):
      # For every band q
      for q in range(self.omega_nq.shape[1]):
        omega = self.omega_nq[n][q]
        # Small frequencies replaced by small numbers
        if omega <= 10**-4:
          omega = 0.001
        H_T = H_T + \
              ( (omega/2.) + self.T*(np.log(1-\
                 np.exp(-omega/self.T))) )
    return H_T/self.omega_nq.shape[1]   


  def get_Cv_T(self):
    """
    Computes the contant-volume Specific heat
    @output
      Cv_T= \del E/\del k_bT
      Note that k_b is already absorbed in units of T.
      You may need to multiply k_b to get proper unit
      NOTE: Cv_T/nkp is returned
    """
    Cv_T = 0.0
    # For every band n
    for n in range(self.omega_nq.shape[0]):
      # For every band q
      for q in range(self.omega_nq.shape[1]):
        omega = self.omega_nq[n][q]
        befac = self.befactor(omega)
        Cv_T = Cv_T + \
              ( (np.square(omega/self.T))*\
                ((1+ (1./befac))*np.square(befac)) )
    return Cv_T/self.omega_nq.shape[1]


  def get_S_T(self):
    """
    Computes the total entropy at T
    """
    S_T = 0.0
    # For every band n
    for n in range(self.omega_nq.shape[0]):
      # For every band q
      for q in range(self.omega_nq.shape[1]):
        omega = self.omega_nq[n][q]
        # Small frequencies replaced by small numbers
        if omega <= 10**-4:
          omega = 0.001
        fac = omega/(self.T*2)
        S_T = S_T + \
              ( omega/(2*self.T*np.tanh(fac)) -\
                np.log(2*np.sinh(fac)) )
    return S_T/self.omega_nq.shape[1]
  

  def get_TS(self):
    """
    Computes the TS so that unit conversion is easy
    """
    TS = 0.0
    # For every band n
    for n in range(self.omega_nq.shape[0]):
      # For every band q
      for q in range(self.omega_nq.shape[1]):
        omega = self.omega_nq[n][q]
        # Small frequencies replaced by small numbers
        if omega <= 10**-4:
          omega = 0.001
        fac = omega/(2*self.T)
        TS = TS + \
              ( omega/(2*np.tanh(fac)) -\
                self.T*np.log(2*np.sinh(fac)) )
    return TS/self.omega_nq.shape[1]
