# Finite Temperature Properties

Everything we've done so far has in fact been for systems at effectively zero temperature. And even at zero temperature, we haven't taken into account the zero-point motion of the atoms and have treated them as purely classical objects residing at fixed positions.

In actuality, there will be an effect from the zero-point motion of the atoms, and as temperature increases the atoms will vibrate with larger amplitude and this can lead to changes in many properties of a material with temperature. You will learn how the thermodynamic properties of a material can be computed from first principles and be able to predict properties such as the total energy and the heat capacity.

Our approach will be to use the type of density functional theory (DFT) and density functional perturbation theory (DFPT) calculations you've already seen, and spend more time analysing the output to produce materials properties. To get the finite temperature properties, we need to get a list of the phonon modes available to the system. The more phonon modes we include in this list, the more accurate the calculation will be, but this comes at the cost of increased computational time. Once the phonon mode list has been obtained, we use it to compute the phonon density of states, which can then be used to obtain the thermodynamical quantity of interest. This will be done using [python](https://www.python.org/).

## Phonon density of states
The phonon density of states is obtained by counting how many phonon modes exist at each frequency. In order to do this, we first need to gather as much information as possible about the phonon modes by calculating the frequency of each mode in a grid of `q` points.


### Step 1. Phonon calculations on a fine-grid

In [Lab08](../lab08/readme.md) you already learned most of the steps you need to follow to perform phonon calculations on a grid. The only difference is the input file for `matdyn.x`.

!!! example "Task 1 - Calculate the dynamical matrix on a fine grid"
    - Copy the lab09 file into your `data` directory and go to the `01_CarbonDiamond` directory. 
    - Using what you learned from [Lab08](../lab08/readme.md), follow the necessary steps to compute the force constant matrix.
    - Take a look at the `05_CD_matdyn-fine.in` input file. You will notice a few differences relative to the previous lab.
    ```bash
    &input
        asr='simple'
        flfrc='CD444.fc'
        flfrq='CD-fine.freq' 
        fldos='CD.dos' 
        nk1=20,nk2=20,nk3=20 #(1)!
        nosym=.true. #(2)!
        dos=.true. #(2)!
    /
    ```
          1. Request a 20x20x20 grid
          2. Do not use symmetry operations to speed up the calculation. This makes it easier to parse the output file
          3. Calculate the density of phonon states

    - Run `matdyn.x` with this input file. It'll take a bit longer than the band calculation as it is explicitly computing without invoking the symmetry.

After `matdyn.x` finishes executing, it will generate two important files: `CD-fine.freq` and `matdyn.dos`. Take a look at the contents of each file. 


!!! example "CD-fine.freq file"
    - The `CD-fine.freq` file contains the phonon frequencies for all the modes at each `q` and is organized as follows
    ```python
    &plot nbnd=   6, nks=8000 /   #(1)!
                0.000000  0.000000  0.000000 #(2)! 
      -0.0000    0.0000    0.0000 1390.4720 1390.4720 1390.4720 #(3)!
              -0.050000  0.050000 -0.050000
      96.9791   96.9791  162.3473 1387.9266 1387.9266 1391.4163
              -0.100000  0.100000 -0.100000
      190.7327  190.7327  321.7228 1380.5533 1380.5533 1393.6345
              -0.150000  0.150000 -0.150000
      278.1010  278.1010  475.2404 1369.1811 1369.1811 1395.4974
        ... 
    ```
          1. The first line tells you there are six bands per `q` point and `8000` (20x20x20) `q` points.
          2. Specifies the `q=(0,0,0)` point.
          3. There are six frequencies on this line, one for each band of the `q=(0,0,0)` point.

Now that you have obtained the phonon vibrational frequencies for a set of `q` points, the phonon density of states is calculated by counting how many phonon modes exist for a given frequency range. For example, if you want to know how many phonon modes exist between frequency $\omega$ and $\omega+d\omega$, this quantity is given by $\rho(\omega)d\omega$.

$$ \text{number of phonon modes between } \omega \text{ and } \omega+d\omega=\rho(\omega)d\omega $$

The density of phonon states is calculated via the following expression:

$$ \rho(\omega) = \sum_{\mathbf{q}\nu}\delta\left(\omega-\omega(\mathbf{q}\nu)\right)$$

where $\delta(...)$ is the Dirac delta function, $\mathbf{q} \nu$ denotes phonon mode $\nu$ at vector $\mathbf{q}$ and $\omega(\mathbf{q}\nu)$ is the frequency of that phonon mode. You can find the result of this calculation in the `CD.dos` file computed by `matdyn.x`.

!!! example "CD.dos file"
    - The `CD.dos` file is organized as follows
    ```python
     # Frequency[cm^-1] DOS PDOS
     -1.3637193144E-05  0.0000000000E+00  0.0000E+00  0.0000E+00 #(1)!
      9.9998636281E-01  1.8451769729E-09  9.2259E-10  9.2259E-10
      1.9999863628E+00  7.3808263018E-09  3.6904E-09  3.6904E-09
      2.9999863628E+00  1.6606947988E-08  8.3035E-09  8.3035E-09
      3.9999863628E+00  2.9523542031E-08  1.4762E-08  1.4762E-08
        ... 
    ```
          1. The first column provides the frequency of oscillation in cm^-1 and the second provides the density of states as states/cm^-1.



We will utilize the density of states to compute several thermodynamic properties using python in the next step.

## Example: Total energy at temperature T
One of the simplest and most important quantities that can be obtained from information about vibrations is the phonon energy. The total energy due to phonons can be written as

$$ E(T) = \int_{-\infty}^{+\infty} E \rho(E)\left(\frac{1}{2} + n_{BE}(E,T)\right)dE $$

where the factor $n_{BE}(E,T) = \frac{1}{\exp(E/k_\mathrm{B} T)-1}$ is recognized as the Bose-Einstein occupation factor.

Notice how the total energy requires considering the contribution of every energy $E$ weighed by the phonon density of states and the Bose-Einstein occupation factor.

### Step 2. Performing the integral using python
This integral will be calculated numerically using python

!!! example "python script"
    - The following script reads in the temperature and the phonon density of states and prints out the total energy at that temperature. Some constants are missing.
    ```python
    import sys
    import numpy as np

    filename = sys.argv[1] 
    temperature = float(sys.argv[2])

    data = np.loadtxt(filename, skiprows=1) #(5)!
    freqs = data[:,0]
    dos = data[:,1]
    constant1 = ???  #(1)!
    constant2 = ???  #(2)!
    dw = freqs[1] - freqs[0] #(3)! 

    total_energy = 0
    for w, rho_e in zip(freqs, dos):

      if abs(w) < 1e-5:
        continue #(4)!

      bose = 1.0/(np.exp(constant1*w/temperature) - 1.0)
      total_energy += constant2*w*(0.5 + bose)*dw

    print(total_energy)
    ```
        1. Replace ??? by the relevant physical constant (or constants) with your choice of units
        2. Replace ??? by the relevant physical constant (or constants) with your choice of units
        3. Numerical integration requires specifying the numerical discretization step
        4. If the frequency is too small, it will cause numerical problems because of a division by a very small number inside the Bose-Einstein occupation factor. Since we're summing over a large number of terms, these ignored terms will have a residual contribution and can be safely ignored
        5. The first row of the `CD.dos` file is a header and should be ignored


Note that this program ignores very small frequencies due to the possibility of dividing by zero. You should fill in the correct numeric value for the constants in this program.


!!! example "Task"
    1. Go back to the lecture notes and fill in the constant missing from the program above
    2. Save this program to a file called `energy.py`
    3. Run this program using the phonon density of states obtained in the previous section, for a temperature of 1K.
    4. Repeat the previous step for a list of ten temperatures from 1K to 100K and plot it using your favourite tools

### Analysing the energy

!!! example "Task"
    1. As the temperature decreases, the energy drops, but eventually it starts stabilizing and never drops below a certain value. Why is this? 
    Answer: As the temperature drops, the phonons lose their thermal energy, but there's still a contribution to their energy coming from the zero-point energy ($\hbar \omega /2$) of each mode 
    2. As the temperature increases, the energy increases linearly with temperature. Why? 
    Answer: For sufficiently large temperatures, the system starts behaving classically, so the amount of energy it can store is proportional to the temperature and the constant of proportionality is the heat capacity.



## Heat capacity at constant volume
Using this type of approach, we are not limited to only computing the total energy. Some key quantities we can calculate are the heat capacity, the entropy and the Helmholtz free energy (For reference, see [Wikipedia](https://en.wikipedia.org/wiki/Quasi-harmonic_approximation) and the reference therein).

The specific heat at constant volume tells us how much energy is required to increase the temperature of a system. It can be obtained from the total energy calculations by differentiating with respect to temperature:

$$C_V(T) = \left(\frac{\partial E}{\partial T} \right )_{const\ V}
 = \int_{-\infty}^{+\infty}dE \left(\frac{E}{k_B T}\right)^2 \rho(E)\frac{\exp(E/k_B T)}{\left[ \exp(E/k_B T)-1\right]^2}$$

!!! example "Task: specific heat"
    - Create a new program `specific_heat.py` that calculates the specific heat of Carbon Diamond by doing the necessary modifications to the `energy.py` program.
    - Calculate the heat capacity for the same temperatures as you did for the energy in the previous task.

### Analysing the heat capacity
At high temperatures, the heat capacity approaches a constant. This is the [Dulong-Petit](https://en.wikipedia.org/wiki/Dulong%E2%80%93Petit_law) law, followed (approximately) by most materials. At low temperatures, the heat capacity is proportional to $T^3$, as predicted by the [Debye model](https://en.wikipedia.org/wiki/Debye_model).

# Summary
In this lab you have learned

- How to compute the phonon density of states
- How to use the phonon density of states to obtain the total phonon energy and the specific heat.