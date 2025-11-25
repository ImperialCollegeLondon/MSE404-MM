Normal Modes and Vibrational Frequencies
========================================

------------------------------ tentative actual text ---------------------------------

This week we continue studying the structural properties of molecules and crystals. In particular, we calculate the vibrational properties of molecules and solids. These properties can be computed using the potential energy surface (PES).  


Normal modes of molecules
-------------------------

### Step 1 


### Step 2 


### Step 3 

!!! example "Task 1: normal modes of methane"

## How quantum espsresso calculates the normal modes of molecules
1. DFPT 
2. Note on choosing the $\Gamma$ point for molecules 
3. Symmetry detection 
4. Stricter ecutrho 



Phonons in crystals 
-------------------------

Force constant matrix
-----------------------

Dynamical matrix 
------------------


## Quantities produced by Quantum espresso 
1. matdyn - how to read 


## How QE calculates the dynamical matrix 


Phonon band structure 
---------------------

## How QE calculates the phonon band structure


------------------------------ old content below ------------------------------------

As we saw in the previous lab, we get easy access to the first derivatives of
the total energy with respect to various parameters in a basic self-consistent
calculation of the density of a system via the Hellmann-Feynman theorem.
However, if we want to access second derivatives, so that we can calculate, for
example, vibrational frequencies, we need to do a little more. This is linked 
to the **2n+1 theorem**, which states that the 2n-th or 2n+1-th order 
derivative of energies would require the n-the order derivative of electronic
wave functions.
 

Vibrational Modes of Methane
----------------------------

We'll start by calculating the vibrational modes of a molecule. This
calculation is set up in 2 parts:

1. A self consistent calculation of the density and wavefunction. You might
   notice this is often the first step in many of the calculations we do.
2. A Density Functional Perturbation Theory (DFPT) calculation of the dynamical
matrix and normal modes. As we have a molecule, we only need look at the gamma
point here.

This calculation is set up in the [`01_CH4`](01_CH4) directory.

- First take a look at the [input file for the scf
  calculation](01_CH4/01_CH4_scf.in). This doesn't have any inputs you haven't
  seen before. There are three things you should take note of though:

1. We've specified the gamma point only by asking for a 1x1x1 unshifted grid,
   whereas previously we've usually explicitly asked for the gamma point. This
   latter option uses optimizations that can make the calculation a bit faster
   if we only use the gamma point, however it does not generate output that is
   compatible with the subsequent DFPT calculation of the frequencies, so we
   instead specify the gamma point in the manner done in the file.

2. We've re-oriented the molecule relative to what you've seen previously, and
   the atomic positions have also been optimized. By specifying the atoms in
   this manner, where the components of the hydrogens along each direction
   have the same magnitude but different signs, the code is much better able
   to detect the symmetry of the molecule. This is very important for
   calculations of the vibrations. If the code understands that all the
   hydrogen atoms are equivalent by symmetry, it only needs to calculate the
   derivatives of the energy with respect to one of the hydrogen positions, and
   then it can populate the full dynamical matrix based on the symmetry of the
   system.

3. We've increased the default value of `ecutrho` which is usually four times
   the value of `ecutwfc`. This can help alleviate some issues with acoustic
   mode frequencies at gamma mentioned below.

### _Task_

- Use this input file to run a calculation with `pw.x`, and save the output
  to a file. Check it worked as expected.


Now take a look at the [input file for the DFPT
  calculation](01_CH4/02_CH4_ph.in). This is as follows:

```
phonons of CH4 (gamma only)
 &inputph
  tr2_ph = 1.0d-15
  asr = .true.
 /
0.0 0.0 0.0
```
- This input file is for the `ph.x` code which uses the output of a SCF
  calculation with `pw.x` to do a DFPT calculation of the dynamical matrix.
  Many more input options are available than we've used here. These are
  described in the `INPUT_PH.txt` file of the quantum espresso documentation.
  Alternatively, you can look up
  [Phonon](https://www.quantum-espresso.org/Doc/INPUT_PH.html)

    - Here the first line is a comment. This can be whatever description you
      like of your input file
    - This is followed by a single section: `INPUTPH`. We've specified two
      variables here.
        - `tr2_ph` which defines the threshold for self-consistency, with a
          default of 1.0d-12.
        - `asr` which turns on the use of the acoustic sum rule which sets the
          three translational modes at gamma to zero (or at least near zero).
          If you take a crystal (or any material) and translate the whole
          system by some amount in real space, the total energy of the system
          remains unchanged. If you have a material in three dimensions, you
          should have three directions where similar translations are possible.
          Therefore, the phonon polarization vectors and frequencies should
          reflect such translational invariance. A material in three dimensions
          should have three phonon modes with zero energy at gamma (as you can
          translate the material along the x, y, and z directions and it does not cost
	  any energy). These phonon modes are referred to as acoustic phonon
          modes. In practice, however, the phonon frequencies are never zero-energy as we
	  use a Fast-Fourier-Transform (FFT) grid in our PW calculations. Other
	  factors, such as the tr2_ph parameter can contribute to acoustic
	  phonon frequencies being non-zero. A simple solution to obtain
	  acoustic phonon modes with zero energy is to impose the acoustic sum
	  rule. The sum rule can be imposed on the gamma point dynamical matrix
	  or on the Force-constant matrix (Force constants are obtained using a
	  Fourier Transform of the Dynamical matrix).  For uniform
	  displacements of all the atoms along x/y/z the net forces acting on
          the atoms are set to zero.
 
    - Following this section, and depending on what input variables have been
      defined, we can give a list of wavevectors at which we want to calculate
      the vibrational frequencies. As with the electronic problem, we are
      modelling a molecule, so only the gamma point makes sense.

### _Task_

- Use `ph.x` to run the calculation using the input file `02_CH4_ph.in` and save
  the output to a file.
- You'll see along with the output file, several files and folders have been
  generated:
    - `matdyn` (the name of this file can be changed in the input file) - this
      stores the dynamical matrix, which can then be passed to several other
      codes for further analysis. Take a look at it now.
        - The top of the file has the atoms, their masses (in internal units)
          and positions.
        - This is followed by the dynamical matrix (3x3 complex numbers for
          each possible pair of atoms - all the imaginary parts are zero at
          gamma).
        - Finally the frequencies of each calculated vibrational mode, both in
          THz and cm-1 are given along with the corresponding wavevector. A
          molecule with 5 atoms will have 15 possible vibrational modes. This
          is how you count vibrational modes/phonons in your system. The 
          same is valid for crystals at any momentum. If you have N atoms 
          within the unit-cell of a crystal, you expect d*N phonon modes
          for every momentum. In all our simulations, d is 3. In a conventional
          material, 3 of these vibrational modes are acoustic i.e. 
          the energy cost associated with the gamma point is zero. On the other
          hand, the rest of the modes (3N-3 of them) cost energy at the gamma point.
          These modes are referred to as optical modes. These modes denote
          relative displacements of the atoms within the unit cell. Both
          acoustic and optical phonon modes can be probed in experiments.
          For example, Raman spectroscopy probes the optical modes and
          is routinely used in studying materials.  

    - `_ph0` -  contains intermediate and restart data for the `ph.x`
      calculation.

- Now let's look through the output of the `ph.x` calculation.
    - Near the top you'll see an analysis of the number of calculations that 
      need to be done given the symmetry of the system. You'll note that first
      3 of the 15 total representations (see below) are computed using the acoustic sum
      rule, the next 3 modes are explicitly computed.  This is why it's so
      important to have the code detect the symmetry of your system correctly.
      Potentially the calculation could have been 5 times slower. Beyond the
      undeniable computational speed up, these representations are quite useful in
      analyzing the phonon modes, such as during phase transitions in materials, and
      Raman activity. They are derived from the so-called **theory of irreducible
      representation of symmetry groups**. Each phonon polarization vector/eigenvector
      denotes a pattern on the crystal lattice. Any arbitrary mode will not obey the
      full symmetry of the crystal (A symmetry of a crystal is an operation that
      leaves the crystal unchanged). The purpose is to group different vibrational
      modes associated with their irreducible representations. One often performs all
      the symmetry operations of the crystal on the phonon polarization vectors and
      categorizes the modes based on the changes caused by the symmetry operations.
      Discussion beyond this is out of the scope of this course.  

    - This is followed by a self consistent calculation for each of the
      representations that need to be calculated.
    - Then we have the frequencies obtained by diagonalizing the dynamical
      matrix.
    - Finally we have a symmetry analysis of the resulting modes, giving
      the mode symmetry and whether they should be detectable by Raman or
      infrared spectroscropy or both.
- You'll likely have noticed we have 3 quite negative frequencies. These are
  the rotational modes of the molecule. It's difficult to get these to come
  out to be zero in the DPFT calculation. In principle this could be enforced
  by a suitable constraint in the code, but this is not available in `ph.x`.
  It's possible in some other codes however.
- The subsequent three modes near zero are the three translational modes, and
  these have been enforced to be close to zero by the `asr` setting in the
  input file.

**Note** as the frequency is calculated using the square root of the
curvature of the energy, a negative frequency is in fact a convention used
to indicate that the frequency is imaginary, and the energy curvature is
negative. This can indicate some instability in the system, so it is good
to ensure you are starting from an optimized structure as you learned how
to do last week.

### _Task_

- Try running the calculation with a higher energy cut-off to see how well
  converged the vibrational frequencies are. Are the vibrational frequencies
  more or less sensitive to cut-off energy than the total energy?

Vibrational Modes of CFC
------------------------

Trichlorofluoromethane, also known as CFC-11 or Freon-11, was formerly used as
a refrigerant. It has accumulated in the atmosphere and is a minor contributor
to the greenhouse effect. This is because, like methane it absorbs radiation
in the thermal IR range where the Earth emits.

### _Task_

- Run the input files in the folder `02_CFC` using pw.x and ph.x, following the same
  process as for methane. Since we need new pseudopotentials (for Cl and F), you will
  need to copy them from `/opt/Courses/MSE404/pseudo` to your own `pseudo` folder.
- CFC-11 has a similar (quasi-tetrahedral) structure to methane, however the ph.x
  calculation takes a lot longer. Why is this? Hint: think about symmetry.
- How do the vibrational frequencies compare to methane? Are they consistent with
  the fact that CFC-11 is a greenhouse gas?
- The cut-off energy is set to be relatively small to speed up the calculation, so 
  that the final results are not well converged. (You may also notice warnings in
  the output file about negative density (rho)). It would take too much time to 
  systematically converge the vibrational frequencies with respect to cut-off, but
  you could try increasing the cut-off slightly. How does the change in vibrational
  frequencies compare to methane?


Phonon Bandstructure of Carbon Diamond
--------------------------------------

Now that you've seen how to calculate the frequencies at a single wavevector,
the next step is to calculate the phonon band dispersion of a crystal. You will
learn how to generalize the theory of phonons in molecules to crystals in one
of the homework.  The directory [`03_CarbonDiamond`](03_CarbonDiamond) contains
a set of input files for carbon diamond.

This calculation has five steps:

1. Perform self-consistent calculation of the density and wavefunction.
2. Calculate the dynamical matrix on a set of wavevectors. We'll call these
   **q-points** from here on, and use **k-points** to refer to electronic
   wavevectors.
3. Fourier transform our grid of dynamical matrices to obtain a set of
   real space force constants.
4. Perform an inverse Fourier transform of the real space force constants to
   obtain the dynamical matrix at an arbitrary q-point. This allows us
   calculate mode frequencies for a fairly dense set of points along lines
   between high-symmetry points in the same manner as for the electronic
   band structure.
5. Generate the plot.

- First examine the scf input file in
  [`01_CD_scf.in`](03_CarbonDiamond/01_CD_scf.in). This is similar to those
  we've seen before. This time we've specified `prefix`. This can be useful if
  you've got different calculations in the same directory. Again, we've
  increased `ecutrho` from its default. The rest is as before.

### _Task_

- Run `pw.x` using this input file and save the output.
    - The other output files generated will be in `CD.save` since we set the
      prefix.
    - Check the output to make sure it completed all right.


- Now examine the `ph.x` input file in
  [`02_CD_ph.in`](03_CarbonDiamond/02_CD_ph.in).
  - This is structured as before, but now:
     - We have specified the same prefix as in the scf input file
     - We've set `ldisp = .true.` which says we're going to calculate a
       grid of q-points.
     - `nq1`, `nq2` and `nq3` define our q-point grid then.

### _Task_

- Run `ph.x` using this input file and save the output. This will take maybe 5
  or 6 minutes to complete, so feel free to read ahead while it's running.
- Once complete you'll see several files have been generated in the
  calculation directory.
    - `matdyn0` has the q-point grid used followed by a list of the q-points
      at which the dynamical matrix has been calculated.
    - `matdyn1`, `matdyn2`, etc is the dynamical matrix at each calculated
      q-point as before.
- Now take a look through the main output file.
    - You'll see it's quite similar to the methane case, but now contains
      a section for each calculated q-point where it figured how how
      many represenations need to be calculated based on the symmetry,
      an scf cycle for each, and the frequencies and mode symmetries
      obtained by diagonalizing the dynamical matrix at that q-point.
    - It's good to look through this file and ensure you're getting
      reasonable numbers for the frequencies before proceeding.


- Next we'll be generating the real space force constants from our grid of
  dynamical matrices using the `q2r.x` code.
- `q2r.x` doesn't have a help file like the other codes. You need to inspect
  the source file in `PHonon/PH/q2r.f90` if you download a copy of the
  quantum espresso source code. The inputs are all described in a comment
  at the top of this file. On the mt-student server you can see this file at
  `/opt/build/quantum-espresso/q-e-qe-6.3/PHonon/PH/q2r.f90`.
- Take a look at the `q2r.x` input file
  [`03_CD_q2r.in](03_CarbonDiamond/03_CD_q2r.in). The contents are as follows:
```
 &input
   fildyn = 'matdyn'
   zasr = 'simple'
   flfrc = 'CD444.fc'
 /
```
- Here we've specified:
    - `fildyn`, which is the name of the dynamical matrix files from `ph.x`,
      where we left it at the default value of `matdyn`.
    - `zasr`, which is the approach used to enforce the acoustic sum rule and
       make the acoustic modes go to zero at the gamma point.
    - `filefrc` is the name of the file in which to output the force constants
       in real space.

### _Task_

- Run `q2r.x` with this input file now and save the output. This will run
  almost instantly.
- The output isn't particularly interesting. The `CD444.fc` file has the
  force constants for each pair of atoms in each direction for each of the
  4x4x4 unit cells in real space.


- Finally we want to use this to generate our normal mode dispersion. We'll be
  doing this with the `matdyn.x` code. As with `q2r.x` this doesn't have a doc
  file describing its input variables. But you can check the comments at the
  top of its source file to get their details. On the mt-student server this
  is at `/opt/build/quantum-espresso/q-e-qe-6.3/PHonon/PH/matdyn.f90`.
- Take a look at our input file
  [`04_CD_matdyn-bands.in`](03_CarbonDiamond/04_CD_matdyn-bands.in). The
  content is as follows:
```
 &input
    asr = 'simple'
    flfrc = 'CD444.fc'
    flfrq = 'CD-bands.freq'
    dos=.false.
    q_in_band_form=.true.
 /
8
 0.000 0.000 0.000 30
 0.375 0.375 0.750 10
 0.500 0.500 1.000 30
 1.000 1.000 1.000 30
 0.500 0.500 0.500 30
 0.000 0.500 0.500 30
 0.250 0.500 0.750 30
 0.500 0.500 0.500 0
```
- Here we're setting:
    - `asr` to again tell it to try to force the acoustic modes at gamma to go
      to zero.
    - `flfrc` to  give it the name of the file with the real space force
      constants from the `q2r.x` calculation.
    - `flfrq` to give it the name of the output file to store its calculated
      frequencies.
    - `dos=.false.` to tell it we're not calculating a density of states
    - `q_in_band_form=.true.` to tell it we want to calculate bands between high
      symmetry points.
    - We finally give it a number and list of high symmetry points with the
      number of points to calculate along each line, in the same way as we did
      for the electronic band structure.

### _Task_

- Run `matdyn.x` using this input file now and save the output. Again this
  is very fast.
- There's very little actual output from the code itself, but it will generate
  the files `CD-bands.freq` and `CD-bands.freq.gp`. Both of which contain
  the frequencies along the lines we requested.
- Finally we want to generate a plot of these frequencies. We could do that
  directly with the previous output: `CD-bands.freq.gp` is a multicolumn
  space separated list of the frequencies in cm-1.
    - It can be easier to use the `plotband.x` tool to help generate a plot.
    - Call this with `plotband.x CD-bands.freq`. This will then ask you for
      an Emin and Emax value - you should pick values equal to or below and
      above the numbers it suggests. Then it will ask you for an output file
      name. Pick "CD-bands-gpl.dat" here. You can then cancel further running
      of this code with `ctrl+c`. Note it has output the location of the
      high-symmetry points along its x-axis.
- Once you've done this, a gnuplot script `plotbands.gplt` has been provided
  that you can use to generate the band structure plot. This is very similar
  to the one used for the electronic band structure, but the location of
  the high-symmetry points along the x-axis has changed as has the name of
  the data file we're plotting. Run this with `gnuplot plotbands.gplt`.

------------------------------------------------------------------------------

Summary
-------

In this lab we have seen

- For a molecule:
    - How to use the `pw.x` code to calculate a converged density and
      wavefunction and then use these with the `ph.x` code to calculate the
      vibrational modes using DFPT.
- For a crystal:
    - How to use the `pw.x` code to calculate a converged density and
      wavefunction.
    - How to use these with the `ph.x` code to calculate the phonon modes on
      a grid of wavevectors.
    - How to transform these to obtain real-space force constants using the
      `q2r.x` code.
    - How to use the real-space force constants to obtain a phonon mode
      bandstructure along lines between high-symmetry points in the
      Brillouin zone.
    - How to plot the bandstructure with gnuplot in a similar way to how we
      plot the electronic bandstructure in a previous lab.

------------------------------------------------------------------------------


### Extra: _Task_

- Calculate the phonon band structure of silicon.
