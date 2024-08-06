Additional Material: Calculating Useful Properties from Total Energies
======================================================================

Once you can calculate a converged total energy for a given structure, you
already know enough to calculate many useful materials properties.

What can we do with total energies?
-----------------------------------

The DFT total energy on its own is not a very useful number, but changes and
differences in the total energy with respect to some input parameter can be
used to calculate many materials properties.

Predicting volume
------------------

If we calculate the total energy of say diamond as a function of volume, the
volume at which the total energy is minimized will be the theoretically
predicted volume.

Let's do this in the directory [`01_diamond_volume`](01_diamond_volume). We
could set up a series of calculations manually, where we generate a set of
input files each with a slightly different value of the input lattice length.
Then we would run each of these calculations, and finally gather the total
energy vs cell volume in a data file. Or, we could write a script that would
do all this for us.

We already have a template file set up as
[`C_diamond_base.in`](01_diamond_volume/C_diamond_base.in). Take a look at
this now. It's set up such that the value of the `A` variable will be replaced
in the automatically generated input files. Note, we're using a fairly high
energy cut-off here also. **In principle once the calculation is set up, we
will need to explicitly test the convergence of the predicted volume versus
energy cut-off.**

The script will need to be slightly more complex than previously, as we'd
managed to get by with just producing integer values in our scripts
previously. Now we'll need to produce floating point values for `A` to get
the resolution we would like. To do this, we can use the
[`bc`](../../misc/linuxcommands/readme.md#bc) command to perform the calculation
of the lattice length for a given input file.

Also Quantum Espresso outputs the unit cell volume (in Bohr cubed) so we can
read this rather than calculating it ourselves.

This gives us the following script:

```bash
#!/bin/bash

template="C_diamond_base.in"
repstr="xxxx"

# The "for" construction we've used previously only handles integers.
# So we set an initial value for A and its delta value as variables.
a1=3.40
da=0.02

for i in {00..10..1}
do
  inp="C_diamond_${i}.in"
  # bc is a calculator that can calculate expressions in the shell. We
  # pass it our calculation with echo and save the result as a variable.
  a=$(echo "$a1 + $i * $da" | bc)
  sed "s/$repstr/$a/" $template > $inp
  pw.x < $inp &> ${inp%.*}.out
done

awk '/unit-cell volume/{vol=$4}
     /^!.*total/{print vol, $5}' *out > etot_v_vol.dat
```

Save this script in the calculation directory and use it to obtain a
file with the total energy versus volume of carbon diamond.

Let's generate a plot of this with `gnuplot` and make sure it looks sensible.
A simple `gnuplot` script [`etot_v_vol.gpl`](01_diamond_volume/etot_v_vol.gpl)
is already in the calculation directory. The contents are as follows:

```
set title "Carbon Diamond, ecutwfc=60 Ry"
set xlabel "Unit Cell Volume (Bohr^3)"
set ylabel "Total Energy (Ry)"
set term pngcairo
set output "etot_v_vol.png"
plot "etot_v_vol.dat"
```
This sets a title and axes labels, sets the output type to png, and sets and
output filename, then plots the data as we've seen before. All these commands
could be entered directly in gnuplot, but it can be easier to save them as a
script if you want to come back in the future and make minor modifications. If
your run this script with `gnuplot etot_v_vol.gpl` you'll see a png file of
the plot has been produced in the directory. You can quickly view this with
`display etot_v_vol.png`.

### _Task_

- What is the calculated volume with the lowest total energy?
- What lattice length does this correspond to?
- Modify the script to calculate 21 points from A=3.40 to A=3.60 inclusive.
- Repeat this last calculation for energy cut-offs of 40 Ry and 80 Ry. Note,
  you should set a different filename for the data file each time or rename
  the previous data file before running this so you can compare them all
  directly..
- Plot all the data for all three energy cut-offs together in gnuplot. Note
  you can plot several data files in gnuplot by separating them with a comma -
  e.g. `plot "data1.dat", "data2.dat"`. Save the output as a png file.

Bulk Modulus of diamond
-----------------------

You may have realised that we could use the results of our previous
calculations to predict the bulk modulus of diamond. The bulk modulus is
proportional to the second derivative of the energy with respect to volume at
the equilibrium volume. Say we approximate our crystal as a harmonic solid,
then we could fit a harmonic expression to it to obtain a value for the
bulk modulus, _K_ in `gnuplot`. We can do that, and plot the fit together
with our data with the following `gnuplot` script.:

```
# Define a function for a simple harmonic equation of state
E(x) = E0 + K*(x - V0)**2/(2*V0)
# It's better to give initial guesses for a good fit. You can make a good
# guess at E0 and V0 from the data plot you generated earlier.
E0=-22.782
V0=74

# The actual fit command. The parameters and errors will be output to your
# terminal and to the file 'fit.log'.
fit E(x) "etot_v_vol.dat" via E0, V0, K

# Let's also produce a plot of the fit along with our data. We can use the
# same settings as before.
set title "Carbon Diamond, ecutwfc=60 Ry"
set xlabel "Unit Cell Volume (Bohr^3)"
set ylabel "Total Energy (Ry)"
set term pngcairo
set output "etot_v_vol_fit.png"

# The fit parameters are set a their optimized values at the end of the fit
# command above, so we can simply plot E(x).
plot "etot_v_vol.dat", E(x)
```

Looking at the fit parameters, we can see that we obtain a value for _K_ and
the equilibrium volume _V0_ directly. The value of _V0_ will be in whatever
units were used for input volumes, so Bohr^3 in our case. The value of _K_
will be in units of Energy/Volume. If we want to compare to the experimental
value we'll need to convert to a more standard unit. 1 Ry is 2.179872325E-18
J, and 1 Bohr^3 is (5.2917721067E-11)^3 m^3. If you don't want to be copy and
pasting numbers to a spreadsheet or similar, another useful way to do
calculations in the terminal is using python. Type `python` to start an
interactive session. Type whatever expressions you want, and type `ctrl+d` or
`quit()` to exit python once you're done.

```
python
>>> K * 2.179872325E-18 / (5.2917721067E-11)**3 # Enter your K value here.
>>> _ * 1e-9 # Convert the last result ('_') to GPa
```

The usual measured value for the bulk modulus of diamond is around 540 GPa.
How close to this are we?

### _Task_

- Convert the fit error to GPa also so we have an estimate of our accuracy.
- Repeat the above fit for the data produced with different energy cut-off
  earlier. How do the values and errors compare?
- A more advanced expression for the equation of state in a solid is the
  [Birch Murnaghan Equation of
  State](https://en.wikipedia.org/wiki/Birch%E2%80%93Murnaghan_equation_of_state).
  See if you can fit this expression to your data. The E(V) expression could
  be written in gnuplot as
  ```
    E(x) = E0 + 9.0/16.0*V0*K*( Kp*((V0/x)**(2.0/3.0) - 1)**3 + ((V0/x)**(2.0/3.0) - 1)**2 * (6.0 - 4.0*(V0/x)**(2.0/3.0)))
  ```
    - Note the expression also depends on the parameter _Kp_ so you'll need
      to add this to the list of fit parameters in the gnuplot `fit` command.
    - How does the error in the fit value of `K` change with this new fit?
    - Generate a plot showing both fits for your most well converged set of
      data along with the data points.

H2 Bond Length
--------------

For simple molecules we can calculate the total energy versus bond length in
the same way. This calculation is set up in the directory
[`02_H2_bond`](02_H2_bond). Again, we're going to use a short script to modify
a single value in a template input file and run a series of calculations. Take
a look at the template input file first. We're not using any new inputs here.
We've placed one atom at the origin, and the second will be moved along the
x-axis.

We'll make a few modifications to the script this time though so that we parse
the data from the output files as they are generated and use the bond-length
value from the script. Save the following to a script in the calculation
directory:

```bash
#!/bin/bash

template="H2_base.in"
repstr="xxxx"

# Again we set an initial value for the x-coordinate of on H-atom and a delta
# value as variables.
hx1=0.700
dhx=0.005

# Empty the file, since we'll be appending to it in the calculation loop.
> etot_v_bl.dat

for i in {00..40..1}
do
  inp="H2_${i}.in"
  # We save the output filename to a variable also.
  out="${inp%.*}.out"

  # Again we use bc to get the atomic position for each input.
  hx=$(echo "$hx1 + $i * $dhx" | bc)
  sed "s/$repstr/$hx/" $template > $inp
  pw.x < $inp &> $out

  # awk is inside the loop this time, and we are appending to the data file
  # after each calculation completes.
  awk -v bl=$hx '/^!.*total/{print bl, $5}' $out >> etot_v_bl.dat
done
```

Take a look at the data file. What bond length minimizes the energy?

Plot this data in `gnuplot`. You may want to copy and modify the script we
used earlier for diamond.

H2 Vibration Frequency
----------------------

As you may have realised, we can use what we already calculated here to find
the H2 vibration frequency. If we assume the variation in energy with bond
length is approximately quadratic (how good an assumption do you think this
is?) we can fit a function to the data in `gnuplot` and find the vibrational
frequency.

We can modify the gnuplot script we used earlier to fit and plot this function
with our data as follows:

```
E(x) = E0 + k*(x-x0)**2
fit E(x) "etot_v_bl.dat" via E0, x0, k

set title "Hydrogen Molecule"
set xlabel "Bond Length (Angstrom)"
set ylabel "Total Energy (Ry)"
set term pngcairo
set output "etot_v_bl_fit.png"

plot "etot_v_bl.dat", E(x)
```

Use this script to make a png plot and inspect it. How good a fit do you think
we have?

### _Task_

- Modify the calculation script to calculate 40 values over a range of
  0.1 Angstrom near the minimum. Also repeat the fit.
    - How does the fit compare to the previous one?
- The fit has produced a value for an effective spring-constant _k_ for the
  H2 vibration. _k_ is in units of Rydberg per Angstrom^2. Find the frequency
  of vibration of the molecule in THz. Note, since we calculated _k_ by fixing
  one atom and moving the other, you will need to use the reduced mass of
  the pair of atoms.
    - What would the frequency be for a molecule made of two deuterium atoms?
