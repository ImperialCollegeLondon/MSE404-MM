# Finite Temperature Properties

Everything we've done so far has in fact been for systems at effectively zero temperature. And even at zero temperature, we haven't taken into account the zero-point motion of the atoms and have treated them as purely classical objects residing at fixed positions.

In actuality, there will be an effect from the zero-point motion of the atoms, and as temperature increases the atoms will vibrate with larger amplitude and this can lead to changes in many properties of a material with temperature. You will learn how the thermodynamic properties of a material can be computed from first principles and be able to predict properties such as the total energy, heat capacity, Helmholtz free energy and Entropy.

Our approach will be to use the type of density functional theory (DFT) and density functional perturbation theory (DFPT) calculations you've already seen, and spend more time analysing the output to produce materials properties. To get the finite temperature properties, we need to get a list of the phonon modes available to the system. The more phonon modes we include in this list, the more accurate the calculation will be, but this comes at the cost of increased computational time. Once the phonon mode list has been obtained, we sum the individual contribution of each phonon mode to get the thermodynamical quantity of interest. This will be done using [python](https://www.python.org/).


## Example: Total energy for specific temperature in diamond
The total energy due to phonons can be written as

$$E(T) = \sum_{\mathbf{q}\nu}\hbar\omega(\mathbf{q}\nu)\left[\frac{1}{2} + \frac{1}{\exp(\hbar\omega(\mathbf{q}\nu)/k_\mathrm{B} T)-1}\right]$$

where the factor $n_{BE}(T) = \frac{1}{\exp(\hbar\omega(\mathbf{q}\nu)/k_\mathrm{B} T)-1}$ is recognized as the Bose-Einstein occupation factor.

Notice how the energy requires summing over every single mode $\nu$ and wave-vector $\mathbf{q}$. The amount of wave-vectors that we need to sum over depends on the choice of $\mathbf{q}$-point grid used to run the DFT calculation. For example, in a $10\times10\times10$ grid, you would need to sum over $1000$ vectors. In general terms, finer grids provide better results, but their computation will also take much longer. To do this calculation, we first need to obtain the list of phonon modes

### Step 1. Phonon calculations on a fine-grid

We've already done most of the work needed to also calculate phonons on a fine-grid in our previous lab ([Lab06/Diamond](../lab06/03_CarbonDiamond)). This can be done following the `q2r.x` calculation by choosing some different input options for `matdyn.x`.  Take a look at the input file [`05_CD_matdyn-fine.in`](CarbonDiamond/05_CD_matdyn-fine.in). The contents are as follows:

```
 &input
    asr='simple'
    flfrc='CD444.fc'
    flfrq='CD-fine.freq'
    nk1=20,nk2=20,nk3=20
    nosym=.true.
    dos=.true.
 /
```

This is quite similar to the band plot, but now we're setting `nosym` to true, choosing a dense q-point grid on which to recalculate our frequencies. Run `matdyn.x` now with this input file. Please do not forget to copy all the things you did in your previous lab on the same material (Copy the files from ([Lab06/Diamond](../lab06/03_CarbonDiamond))
directory into [Lab07/CarbonDiamond](lab07/CarbonDiamond). It'll take a bit longer than the band calculation as it is explicitly computing without invoking the symmetry. After it finishes, it will generate the following file: `CD-fine.freq`. Take a look at the contents of this file. It is organized as such

```
freq1
freq2
freq3
```

We will utilize this file to compute several thermodynamic properties using python in the next step.

### Step 2. Summing over modes in python
To calculate the total energy due to phonons using python, we will create a simple script that reads in the temperature and a file with the list of frequencies and prints out the energy.

```python
import sys
frequencies = sys.argv[1]
temperature = sys.argv[2]
energy = 0
for frequency in frequencies:

  if abs(frequency) < 1e-5:
    continue

  x = frequency/temperature
  bose = 1.0/(exp(x) - 1)
  energy += x*(0.5 + bose)

print(energy)
```

Note that this program ignores very small frequencies due to the possibility of dividing by zero.

Thermodynamic properties
=======================
Some key quantities are (For reference, see
[Wikipedia](https://en.wikipedia.org/wiki/Quasi-harmonic_approximation) and
 and the reference therein):

Bose-Einstein distribution
--------------------------

$n_{BE}(T) = \frac{1}{\exp(\hbar\omega(\mathbf{q}\nu)/k_\mathrm{B} T)-1}$

Total Energy due to phonons
---------------------------
Total energy due to phonons within harmonic approximation can be 
written as,

$$E(T) = \sum_{\mathbf{q}\nu}\hbar\omega(\mathbf{q}\nu)\left[\frac{1}{2} + \frac{1}{\exp(\hbar\omega(\mathbf{q}\nu)/k_\mathrm{B} T)-1}\right]$$

Constant volume heat capacity
-----------------------------
Specific heat at constant volume can be obtained from the total 
energy calculations:

$$C_{V} = \left(\frac{\partial E}{\partial T} \right ) = \sum_{\mathbf{q}\nu} k_\mathrm{B} \left(\frac{\hbar\omega(\mathbf{q}\nu)}{k_\mathrm{B} T} \right)^2 \frac{\exp(\hbar\omega(\mathbf{q}\nu)/k_\mathrm{B} T)}{[\exp(\hbar\omega(\mathbf{q}\nu)/k_\mathrm{B} T)-1]^2}$$

Helmholtz free energy
---------------------
To compute the Helmholtz free energy, we need the partition function,
$Z$.

$$Z = \exp(-\varphi/k_\mathrm{B} T) \prod_{\mathbf{q}\nu} \frac{\exp(-\hbar\omega(\mathbf{q}\nu)/2k_\mathrm{B}T)}{1-\exp(-\hbar\omega(\mathbf{q}\nu)/k_\mathrm{B} T)}$$


$$H(T) = -k_\mathrm{B} T \ln Z = \varphi + \frac{1}{2} \sum_{\mathbf{q}\nu} \hbar\omega(\mathbf{q}\nu) + k_\mathrm{B} T \sum_{\mathbf{q}\nu} \ln \bigl[1 -\exp(-\hbar\omega(\mathbf{q}\nu)/k_\mathrm{B} T) \bigr]$$

Entropy
-------
Entropy, $S$ can also be computed:
$$S = -\frac{\partial H}{\partial T} = \frac{1}{2T} \sum_{\mathbf{q}\nu} \hbar\omega(\mathbf{q}\nu) \coth(\hbar\omega(\mathbf{q}\nu)/2k_\mathrm{B}T)-k_\mathrm{B} \sum_{\mathbf{q}\nu} \ln \left[2\sinh(\hbar\omega(\mathbf{q}\nu)/2k_\mathrm{B}T)\right]$$


Note that the temperature dependence in all these quantities are
determined by the Bose-Einstein distribution.


We have implemented these codes using python. For example, you will
find a folder [`Thermodynamics`](03_CarbonDiamond/Thermodynamics)
containing a file `thermo.py` that has all these quantities you need. It reads
the phonon bands calculations at a fine-grid and can compute several
thermodynamic properties. For the implementation of total energy due to
harmonic phonon, look up the function, *get_E_T* inside the `thermo.py` file:

 ```
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
```

You can compute the total energy at a given temperature by running
the `compute.py`: `python compute.py`

Try understanding how the calculations are performed and how are
they implemented. Note that, all the calculations internaly
converts temperature to meV (i.e. $T \to k_{B}T$) and phonon
energies to meV (from cm$^{-1}$).

### _Task_

- Run the calculations at a temperature 20 K.
- Compute the temperature dependence of $E, C_{V}, H, S$ for 
  several temperatures, ranging from 10 K to 1000 K in steps
  of 20 K. What happens to specific heat at low-temperature?
  You will see a lot more details on your homework.
  **Hint:** Write a *for* loop to do this.
- Plot these data using matplotlib.
- Try increasing the grid-size from $20\times20\times20$ to a larger number and
  try reducing as well. What happens?


**NOTE:**  An important contribution in the total energy,
entropy, etc. are missing in the above calculations: the
contribution without the phonons. For example, the total energy of
a material at a given temperature is, $$ E(T)= E_{ph}(T) + E_{DFT}$$,
where $E_{DFT}$ is the contribution without phonon (sometimes,
referred to as lattice energy). Another important point to note is
that they are all dependent on the volume of the material.
Similarly, the Helmohltz free energy can be written as $H = E_{DFT}+ E_{ph} - TS$. A thermodynamic state is described by two
independent parameters, let's take them as $T$ and $V$. The
free-energy, $H(T,V)$. So, in principle, one should compute the
free-energy for several volumes at any temperature or vice-versa,
to represent the thermodynamic state correctly. This leads us to
something called, **quasi-harmonic approximation**. This
approximation is a harmonic-phonon-based model used
to describe volume-dependent thermal effects, such as
the thermal expansion of a material. This approximation assumes
that the harmonic approximation holds for every value of the
lattice constant, viewed as an adjustable parameter. 

