Finite Temperature Properties
=============================

Everything we've done so far has in fact been for systems at effectively zero
temperature. And even at zero temperature, we haven't taken into account the
zero-point motion of the atoms and have treated them as purely classical
objects residing at fixed positions.

In actuality, there will be an effect from the zero-point motion of the
atoms, and as temperature increases the atoms will vibrate with larger
amplitude and this can lead to changes in many properties of a material with
temperature. You will learn how the thermodynamic properties of a material
can be computed from first principles.

Our approach will be to use the type of density functional theory (DFT) and
density functional perturbation theory (DFPT) calculations you've already seen
, and spend more time analysing the output to produce materials properties.


This will involve some reasonably serious numerical calculations, of the kind
that would be quite difficult to do without using some form of mathematical
software. We will use [python](https://www.python.org/) to compute the
thermodynamic properties.


We'll be mainly relying on the following libraries:

- [numpy](https://www.numpy.org/)
    - This allows us to easily work with arrays of data.
- [scipy](https://www.scipy.org/)
    - This gives us many numerical analysis tools.
- [matplotlib](https://matplotlib.org/)
    - This allows us to easily generate plots.


We expect many of you may not be familiar with python, but if you've used
Matlab or Mathematica, you should find the process here somewhat similar, but
with slightly different syntaxes. Python is a very powerful language that you
are most likely to use in future. We will also introduce you to *functions* and
*classes*. Functions do specific tasks that you wish to perform. For example, 
you provide some input, a function does some calculations and returns you an
output. Therefore, with a *function*, you primarily interact (through input  and
output). On the other hand, a *class* allows you to not only create your own
data type but also interact with it. You will find out that *class* is
significantly more re-usable.


Let's introduce you to a very simple *class*. It has two aspects:
implementation, and interaction/usage. Let's look at implementation
of the position of an atom [position](./class/):
```
# Indentation has to be consistent
# We are using two spaces as indentation.

# Creating a position class
class Position(object):
  # Special method __init__ to initialize your data
  # attributes
  # Note as opposed to normal function, __init__ contains
  # self;
  # self: parameter to refer to an instance of a class
  # x,y: what you provide while creating this class/calling it
  def __init__(self, x, y):
    # self.x or self.y: Look for x/y that belong to this class
    self.x = x
    self.y = y

  # Methods that you can use to compute things
  # Following method computes distance of (x,y) from origin
  def get_dist_from_origin(self):
    # Note how the x/y values are called (with self.)
    dist = (self.x**2.0 + self.y**2.0)**0.5
    return dist

  # You can add as-many methods as you want
```
This is the content of `position.py`. Now, let's see how you can
use it. This class can be called to compute the distance of a
point from origin. To do this, run the code `run.py` using :
`python run.py` . Now, take a look inside `run.py` and try to
understand how the calls are made. 

```
# Import the class Position
from position import Position

pos = Position(2,3)
# Access data from class
print("x=%.6f"%(pos.x))
print("y=%.6f"%(pos.y))
# Compute distance by calling the function from class
distance = pos.get_dist_from_origin()
print("distance from origin(0,0)=",distance)
```

Thermodynamic properties of the material involve understanding the
phonon occupation number (or, Bose-Einstein distribution), total
energy within the harmonic approximation, specific heat at constant
volume, Helmholtz free energy, Entropy, etc. In this lab, you will
try to learn to code some of these properties by building on some
existing code for Diamond.


Phonon calculations on a fine-grid
==================================
As you might have noticed, we've already done most of the work
needed to also calculate phonons on a fine-grid in our previous lab
([Lab06/Diamond](../lab06/03_CarbonDiamond)). This can be done
following the `q2r.x` calculation by choosing some different input
options for `matdyn.x`.  Take a look at the input file
[`05_CD_matdyn-fine.in`](CarbonDiamond/05_CD_matdyn-fine.in). The
contents are as follows:

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

This is quite similar to the band plot, but now we're setting
`nosym` to true, choosing a dense q-point grid on which to
recalculate our frequencies. 

Run `matdyn.x` now with this input file. Please do not forget to
copy all the things you did in your previous lab on the same
material (Copy the files from ([Lab06/Diamond](../lab06/03_CarbonDiamond))
directory into [Lab07/CarbonDiamond](lab07/CarbonDiamond). It'll take a bit
longer than the band calculation as it is explicitly computing without invoking
the symmetry. The file it generates: `CD-fine.freq`. We will utilize this file
to compute several thermodynamic properties using python.


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


Thermal Expansion of a Solid
----------------------------
$$ \alpha = (\frac{1}{V} \frac{\partial V}{\partial T})$$ at a constant
pressure. Note that the volume, $V$ is obtained by minimizing Gibbs
free energy, $G=H+PV$. Assuming a zero-external pressure, we can
minimize H at every temperature and obtain the thermal expansion
coefficient. You will find this in the [Si_TEC](Si_TEC). 

### _Task_
The calculations has three steps:
1. Run the total energy and phonon calculations for several volumes.
2. Compute the Helmohtz free energy for any temperature for these volumes.
3. Obtain the minimum for for each temperature and compute V vs. T.
4. Write a small python code to compute $\alpha$ using the scripts
provided in [TEC](Si_TEC/TEC) folder.

The first-step can be established using running several
calculations`bash EvsV.sh`. After this step, go to
[TEC](Si_TEC/TEC) directory and run `compute.py`.If you recall, it
took around 5 minutes to calculate the phonon band-structure for a
single volume last time. Doing this 6 times would take half an
hour. For this reason, we've included the dos files that would be
generated. **You don't need to run the calculation here and can
proceed to the analysis** [TEC](Si_TEC/TEC)

------------------------------------------------------------------------------

Summary
-------

In this lab you have seen:

- For a solid we used
  - python code to compute several thermodynamic properties computed using DFT+DFPT calculations.
