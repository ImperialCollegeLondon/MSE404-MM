Additional Material: Optical Properties and Time-Dependent Density Functional Theory
====================================================================================

Among the most useful class of properties to be able to predict for a molecule
or crystal are its optical properties. These can be used to find the
frequencies of optical radiation that will be absorbed or emitted, based on
its electronic structure. This can be used to find, for example, what colour
a dye molecule will have in a solvent.

There are several ways these properties could be calculated starting from
density functional theory. The method we'll look at in this lab is
time-dependent density functional theory (TDDFT). A TDDFT code could be used
to calculate the evolution of the electrons under the effect of the
oscillating electric field associated with the presence of a photon. This
would be a fairly intensive calculation, and would need to be repeated at many
different energies to build up the full optical spectrum.

The set of codes that comes with the quantum espresso package, turboTDDFT,
actually use time-dependent density functional *perturbation* theory (TDDFPT)
to calculate optical spectra of molecules. This calculates the response of the
system to the oscillating electric field associated with an incoming photon
using perturbation theory (polarizability), in a similar way to how DFPT finds
the response of the system to perturbations of the atomic positions. It then
uses this to find the susceptibility as a function of photon energy. The code
also uses an approach known as the Liouville-Lanczos method which allows it to
calculate the full optical spectrum for a wide energy range, at a very
moderate cost.

Typically if you want to understand what energy of photon will be absorbed,
you might expect you would need to know the difference between the energy of
the occupied and empty states: when a photon is absorbed it excites an
electron from an occupied state to an empty state higher in energy by an
amount equal to the photon energy. The turboTDDFT codes use a clever approach
to avoid the need for including a large number of additional empty states in
the calculation. This involves instead using the projection operator on to the
occupied states. Details of this are given in a paper outlining the methods
used in the turboTDDFT code. This is distributed with the quantum espresso
package. You can find it on the mt-student server in
`/opt/share/quantum-espresso/doc-6.3/turboTDDFT-CPC.pdf`. There are also
slides on the quantum espresso website outlining the approach and making a
comparison to the more basic TDDFT method at
https://www.quantum-espresso.org/resources/tutorials/shanghai-2013/TDDFT_Talk_Gebauer.handouts.pdf

Despite the various advantages offered by the approaches used in this code,
this will still be the most intensive calculation you'll have done so far.
However it is much faster than would be the case if we needed to manually
look at many photon energies and including many empty states in our DFT
calculation.

You should also keep in mind that this approach is for _fixed nuclear
positions_ and does not include the effect of interatomic vibrations which may
also couple to an incoming photon, either directly if the energy range is
similar, or indirectly if, for example, a photon creates an excited electron
and vibration simultaneously. Generally by looking in a particular energy
range, the dominant effects will be from one mechanism or another. For
example, the homo-lumo gap in methane is around 10 eV in LDA-DFT, while the
highest energy atomic vibrations are around 0.4 eV. This means that if we're
looking at the effect of photons in the 10-50 eV range (vacuum UV to extreme
UV), the effect of the atomic vibrations will be minor. But we don't
immediately have the spectrum across all energy ranges.

As with some of the other more advanced calculations we looked at, we'll be
doing this calculation in several stages:

1. A self consistent calculation of the molecule is performed. We can do this
   in the same way as previously. No special inputs are needed.
2. We'll use the `turbo_lanczos.x` code which does the TDDFPT calculation
   using the aforementioned Liouville-Lanczos method. This is the first, and
   more computationally intensive step in calculating the polarizability. It
   calculates the components of a tridiagonal matrix that will be used in a
   subsequent code, for each optical polarization direction.
3. We'll use the `turbo_spectrum.x` code to take the tridiagonal matrix
   generated in the previous step, generate the polarizability and from it
   find both the susceptibility tensor and the oscillator strength as a
   function of photon energy.
     - The susceptibility tensor relates the electric field to the induced
       electric polarization
     - the oscillator strength gives the probability of absorption or emission
       of photons.

Optical Spectrum of Methane
---------------------------

The directory [`CH4`](CH4) contains a set of example inputs to perform
this calculation for a methane molecule.

- Start by examining the `pw.x` input file
  [`01_CH4_scf.in`](CH4/01_CH4_scf.in). You'll see that we have set this up as
  usual. We have specified both the `prefix = CH4` and `outdir = './out'` as
  this calculation generates many intermediate files, so it will be easier to
  keep them all together. We have also set the positions such that the C atom
  is at the centre of the cell. This will make some subsequent visualization a
  little easier.
- Run `pw.x` with this input file and examine the output to ensure it worked
  as expected.
- Next take a look at the `turbo_lanczos.x` input file
  [`02_CH4_tl.in`](CH4/02_CH4_tl.in). This is fairly short. As usual, we have
  accepted the default values for most of the input parameters. You can see
  the full list of inputs in the `INPUT_LANCZOS.txt` file in the quantum
  espresso documentation directory. We have specified two sections in the
  input file:
    - `LR_INPUT` where we set the `prefix` and `outdir` to match the scf
      calculation, and
    - `LR_CONTROL` where we set
        - `itermax = 400` which tells it how many elements of the Lanczos
          chain to calculate (elements of the tridiagonal matrix). While in
          principle more elements will give us a more accurate result, the
          subsequent code can extrapolate these to a very high number. Usually
          on the order of 500 to 1000 is sufficient to calculate explicitly,
          though we'll look at how this converges later. It is connected to
          how well the spectral features can be resolved in your calculation.
        - `ipol = 4` which tells it to calculate the response in all 3
          directions: x, y and z.
- Run `turbo_lanczos.x` with this input file. It will likely take somewhere
  around five minutes on the mt-student server. If you follow the output as it
  is generated, you will see it calculating each of the 400 requested
  iterations in turn for each of the three polarization directions.
- Finally take a look at the `turbo_spectrum.x` input file
  [`03_CH4_ts.in`](CH4/03_CH4_ts.in). We need to set a few more things here.
  As before, you can see a full description of all the input parameters that
  can be used with this code in `INPUT_SPECTRUM.txt` in the quantum espresso
  documentation directory. It has a single section: `LR_INPUT`. Again we set
  the prefix and outdir to match the previous calculations. Then we have:
    - `itermax0` which needs to match the number of iterations calculated
      explicitly by `turbo_lanczos.x`.
    - `itermax` which says how many iterations in total to calculate. Those
      not explicitly calculated are generated using an extrapolation scheme.
    - `extrapolation` sets the approach used in the extrapolation.
    - `espil` sets a broadening to use (in Ry). Very small wiggles will appear
      in the output spectrum if this is too small, and if it is too large it
      may obscure important features.
    - `units` sets the output units. We use `1` here we means the output
      energies are in eV (rather than Ry by default).
    - `start`, `end` and `increment` specify the range of energies to output
      and the spacing between them (in the units specified by the `units`
      input).
    - `ipol` specifies the polarization direction. We set this to `4` again
      for all three cartesian directions.
- Run `turbo_spectrum.x` with this input file. This is quite a quick
  calculation, and should finish in 10 to 20 seconds. The output file mainly
  contains the stages of the calculation.
- The main output is instead in the file `CH4.plot_chi.dat` which will be
  generated. This lists each energy value and its corresponding susceptibility
  tensor and oscillator strength.

Take a look at the file `CH4.plot_chi.dat` now. You'll see the file is not
really amenable for plotting.

- First we can extract the oscillator strength with the following awk command:

```bash
awk '/S\(E\)=/{print $2, $3}' CH4.plot_chi.dat > ostrength.dat
```
This will print the second and third character on every line containing the
text `S(E)=` (we need to use a `\` to escape the parentheses). This extracts
the oscillator strength as a function of energy (in eV since we set `units =
1` in the `turbo_spectrum.x` input file).

- Plot the file containing the oscillator strength as a function of energy in
  gnuplot so you can see where the peaks in absorption are.

It's also interesting to look at the susceptibility. As the system is
isotropic, the polarizability, and hence susceptibility in each of the three
Cartesian directions is the same. This means looking at the `_11` component
alone is sufficient. This is a complex number, and we should look at both
the real and imaginary parts.

- Parse the energy, the real and the imaginary part of the susceptibility
  (from the 11 component of chi in the `plot_chi.dat` file) and plot both the
  real and imaginary parts together in gnuplot.
    - If you parse this to a file with three columns: x, y1, y2, you can plot
      them together in gnuplot as
      `plot "file" with lines, "file" using 1:3 with lines`.
    - You can plot several plots together by separating them with commas, and
      you can use the keyword `using` to specify the columns to plot (it is
      `1:2` by default).

Charge Density Response
-----------------------

It can also be interesting to examine how the charge density actually varies
in response to a photon at a given energy (and hence frequency) and
polarization. In particular it's useful to look at the response to photons at
the energy where the highest peak in the oscillator strength is.

This can be done using an additional flag and input section in a
`turbo_lanczos.x` calculation, following a previous `turbo_lanczos.x`
calculation for all polarization directions.

Take a look at the input file [`04_CH4_tl.in`](CH4/04_CH4_tl.in). In this
file we have set the `LR_CONTROL` section inputs:

- `itermax = 400` which matches what we had previously.
- `charge_response = 1` which turns on the calculation of the charge response.
- `ipol = 1` which set the polarization of the calculated charge response to
  be along the x-direction.

Then we have added an addition section `LR_POST`, which is needed whenever
`charge_response` is set to 1. Here we have set the following:

- `omeg = 1.08` which set's the energy of the photon for which we want to
  calculate the response. This is in Ry, and corresponds to 14.7 eV which is
  roughly where the peak of the oscillator strength is for this level of
  (under)convergence.
- `epsil = 0.02` this is an energy broadening to use when calculating the
  response (in Ry).
- `w_T_npol = 3` which says that we calculated the response to all three
  polarization directions in our previous calculation.
- `plot_type = 3` which selects the output file to be in the gaussian cube
  format. This can be opened with xcrysden.
    - Note: it seems there may a bug with this output for the version compiled
      on the mt-student server such that for particular choices of energy
      cut-off, box-size and output format, the output file is cut off before
      the end and can't be plot. For the values chosen here, it should work
      hopefully work correctly.

Run `turbo_lanczos.x` with this input file. It'll take a couple of minutes to
finish. Once completed you'll see *either* a pair of files named
`CH4-absorbtive-pol1.cube` and `CH4-dispersive-pol1.cube` have been created,
*or* a single file named `CH4-summed-rho-pol1.cube`. The former are generated
where the code determines the energy `omeg` (including the broadening `epsil`)
you are calculating at corresponds to a resonance where there is a peak in the
spectrum, and the latter is output otherwise. In our case we should be within
the broadening of a resonance peak and so two files should be output.

Open `xcrysden`, which we used to look at crystal and molecular structures in
the extra section at the end of lab 2, (don't forget to load the module). Use
this to open one of the generated files by navigating through File -> Open
Structure -> Gaussian98 Cube File, and selecting the file. You can look at the
regions where the charge density is increased and decreased in response to the
oscillating electric field of the photon by going to Tools -> Data Grid, click
OK. Then select an appropriate Isovalue in the menu that appears, tick the
`Render +/- isovalue` box, and click `Submit`. To find a good isovalue you can
repeat this process with different values within the listed maximum and
minimum on the grid till you see something interesting; the menu will stay
visible after you click submit.

Convergence
-----------

The calculated spectrum is quite sensitive to both the box size, the
plane-wave energy cut off, and the number of elements of the Lanczos chain
explicitly calculated, as set by the `itermax` variable in the
`turbo_lanczos.x` input file. (And subsequently used as `itermax0` in the
`turbo_spectrum.x` input file). These are all somewhat somewhat underconverged
in the input file given to you earlier so that the calculation would complete
in a reasonable time. And are somewhat connected to each other: e.g. the
spectrum may be converged at a higher value of `itermax` for a higher value of
`ecutwfc` so that testing convergence can be time consuming.

Underconverging your spectra can lead to the appearance of spurious peaks and
shoulders and to changes in the locations of different peaks. So you will
always need to explicitly check the convergence of your spectra with respect
to these parameters.

As it would take quite a few heavy calculations to for you to test this, the
oscillator strengths for methane for various box sizes, energy cut-offs and
values of `itermax`
have been pre-calculated in the [`CH4_convergence`](CH4_convergence)
directory. For the various files here, in the filename the box dimension in
Angstrom is indicated by the number following `A`, the plane-wave energy
cut-off in Ry is indicated by the number following `e`, and the value of
itermax is indicated by the number following `i`.
Try plotting and comparing combinations of these in gnuplot to see how the
spectrum changes as these parameters are changed.

> Note: if you are testing the convergence with respect to `itermax`, you can
generate output from `turbo_spectrum.x` where `itermax0` is any value _up to_
the value chosen for `itermax` in the `turbo_lanczos.x` calculation. So for
example, if you have calculated a spectrum with `itermax = 500` in
`turbo_lanczos.x`, you can quickly generate output for say 400 by setting
`itermax0 = 400` in `turbo_spectrum.x`, rather than doing a full recalculation
with a lower value of `itermax` in `turbo_lanczos.x` which would be far more
time consuming.

-------------------------------------------------------------------------------

Summary
-------

In this lab you have seen how to use the turboTDDFT codes for performing and
analysing TDDFT from the quantum espresso package to

- calculate the optical spectrum of a molecule by
    - Performing a self-consistent DFT calculation with `pw.x`.
    - using the `turbo_lanczos.x` code to perform a TDDFPT calculation as
      the first step in calculating the polarizability.
    - using the `turbo_spectrum.x` code to generate the polarizability and
      find both the susceptibility tensor and oscillator strength as a
      function of photon energy.
- You have also used `turbo_spectrum.x` to calculate the charge density
  response of a molecule and visualised it with `xcrysden`.


-------------------------------------------------------------------------------

*Extra - Li2*
-------------

If you want to practice this type of calculation a bit more, the directory
[`Li2`](Li2) contains a pseudopotential for lithium. You can use this
to do the following set of calculations:

- Find the optimal bond length for an Li2 molecule by relaxing the atomic
  positions as shown in previous labs.
- Calculate the oscillator strength as a function of energy.
- You will likely need to play around with the values of the energy range to
  use in the `turbo_spectrum.x` calculation to find the best values to
  show the main peaks.
