Additional Material: The Tetrahedron Method for Densities of States
===================================================================

You saw in [lab 4](../../../lab04/readme.md) how to calculate densities of
states using broadening, in this section we look at how to use the tetrahedron
method.

### Density of states using the tetrahedron method

To use the tetrahedron method this within quantum espresso, you must use
the tetrahedron occupation scheme for the DFT calculation, and then don't
set any `degauss` value in the `dos.x` input file. This is outlined in the
`INPUT_DOS.txt` help file.

An example where we have modified the diamond input files from [lab
4](../../../lab04/readme.md) to calculate the density of states using
tetrahedron integration is given in the [`01_diamond_tet`](01_diamond_tet)
directory. 

If you run the same three steps for calculating densities of states using
broadening again here, you'll once again produce a `pwscf.dos` file. If you 
check the output file from the `dos.x` calculation you'll also see the note
`Tetrahedra used`, whereas the previous output notes details of the broadening.

Additionally if you look at the header of the `pwscf.dos` file you'll see a
different value for the Fermi energy than we obtained in the previous
calculation. This would mainly come from a small underestimation that happens
in the tetrahedron method. You can also see this in the integrated DOS in
comparison to the broadening case. For this reason a separate calculation of
the Fermi energy may often be useful. Note, you'll get more well defined
features and a smoother curve systematically as you increase the k-point
sampling density. The calculation will take a lot longer; if you use e.g. a
60x60x60 grid in the nscf calculation it takes about 15 mins run in serial,
and the DOS calculation takes about 20 minutes, though this may be worth it if
producing an image for a paper or report. We've tried to keep run time for
example calculations provided down to a couple of minutes or less.

### _Task_

- Plot the density of states and compare it that obtained using Gaussian
  broadening in [lab 4](../../../lab04/readme.md).
