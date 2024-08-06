Additional Material: Examining the Charge Density and Wavefunction
==================================================================

In this lab we'll look at various ways of visualizing the output from your DFT
calculations. To do this we'll be using several additional codes from the
Quantum Espresso package.

The Charge Density
------------------

To start with we'll try to visualise the charge density we have calculated for
a methane molecule. Take a look at the directory
[`01_chargedensity/01_methane`](01_chargedensity/01_methane). This is contains
an [input file](01_chargedensity/01_methane/CH4.in) for methane exactly as
we've seen before. Note, we've set `disk_io = 'low'` which is the default
value for an scf calculation (i.e. we could have omitted this and get the same
output), as we want to keep the charge density file for analysis.

- Run `pw.x` with this input file now, and check the output to make sure it
  all worked as expected. You should have both a `pwscf.save` directory
  containing the charge density file and other output, and a `pwscf.wfc`
  wavefunction file.

Now we want to post-process this output into something we can visualise more
easily. To do that we'll be using the `pp.x` code from the Quantum Espresso
suite. This code can read the output files produced by `pw.x`, extract
whatever quantity you're interested in and generate output compatible with
various visualisation programs. There is a help file similar to the `pw.x`
available in the Quantum Espresso documentation folder called `INPUT_PP.txt`.

### 3D visualization with XCrySDen

There is a `pp.x` input file for methane called
[`CH4_pp.in`](01_chargedensity/01_methane/CH4_pp.in), that contains the
following:

```
 &INPUTPP
    filplot = 'methane_charge'
    plot_num = 0
 /

 &PLOT
    filepp(1) = 'methane_charge'
    iflag = 3
    output_format = 5
    fileout = 'CH4.rho.xsf'
 /
```

As you can see, the file format is quite similar to that used by `pw.x`. It
contains two sections:

- `INPUTPP` which has inputs controlling what data will be extracted from the
  `pw.x` output.
    - If we had chosen non-default directories and names for the output files,
      we would need to se them here.
    - `filplot` specifies a filename that will be used to save the extracted
      data.
    - `plot_num` takes an integer option that specifies what quantity will be
      extracted. 0 gives us the electronic density. There are many options
      here, such as potentials and local density of states. See the help
      file `INPUT_PP.txt` for details, including info on additional variables
      that apply depending on this option.
- `PLOT` has inputs controlling how the extracted data will be output for
   visualisation.
   - `filepp(1)` gives the name of the first file to read (it's possible to
     read several files and combine them in various ways, such as to
     generate charge density difference plots).
   - `iflag` gives the dimensionality of the output. We choose 3 here for a
     3D plot.
   - `output_format` gives the type of format you want to have output, usually
     this is determined by what application you'll be using to visualise the
     data. It should also be compatible with the value of `iflag`. We've
     selected 5 which will give us output compatible with `xcrysden`.
   - `fileout` specifies the name of the plot file. `xsf` is a default
     `xcrysden` extension.

We can run `pp.x` in the same way as `pw.x`:

```bash
pp.x < CH4_pp.in &> CH4_pp.out
```

You'll see both a `methane_charge` and a `CH4.rho.xsf` file have been
generated in the calculation directory. The charge file is in a binary format
so we can't tell much about it. The `.xsf` file is in a text format so we
can examine it (and modify it if we want) in a text editor. Take a look and
you'll see the file simply defines the crystal lattice and basis, then has a
section with the datagrid as a mesh of points.

Launch `xcrysden` (recall you'll need to load the xcrysden module first with
`module load xcrysden`), and load the file `CH4.rho.xsf`. You'll see a big
cube, which has the carbon at each corner and whichever bound hydrogen
falls within the cube nearby. This isn't ideal, but let's carry on for the
moment. There's no sign of the charge density yet. To enable this, go to
"Tools" -> "Data Grid". You can click "OK" on the menu that appears, and then
you'll see a menu that controls the data grid plot appearance. Try entering
say "0.1" for the "Isovalue" and click "Submit". If you zoom into a carbon
you'll see the cloud of charge around it. Unfortunately, since we have the
carbon at the origin, this plot is a little hard to see.

### _Task_

- In a new directory `01_methane_shift`, shift the molecule such that the
carbon atom is at the centre of the box in the `pw.x` input file.
- Run the `pw.x` and `pp.x` calculations as before.
- Visualise the output in `xcrysden`.
    - Disable the display of the crystal cell.
    - Plot the data grid at an Isovalue of 0.27. You should also enable
    transparency, and increase the degree of the tricubic cpline to 2 or 3.
    - Under "Modify" reduce the Ball Factor to reduce their size.
- When you're happy with how it looks, save an image by first checking the
  "File" -> "Print Setup" menu. It's usually worth enabling anti-aliasing
  here. Then go to "File" -> "Print" and save the output as a png image.
- Why is there no density around the carbon atom?

### 2D visualization with GnuPlot

The other way to look at this would be as a 2D contour plot, or heat map.
We can get 2D output suitable for gnuplot by setting the input for `pp.x` as
in [`CH4_pp_gp.in`](01_chargedensity/01_methane/CH4_pp_gp.in):
```
 &INPUTPP
 /

 &PLOT
    filepp(1) = 'methane_charge'
    iflag = 2
    output_format = 7
    fileout = 'CH4.rho.gpl'
    e1(1) = 0.4, e1(2) = 0.0, e1(3) = 0.0
    e2(1) = 0.0, e2(2) = 0.0, e2(3) = 0.4
    x0(1) = -0.2, x0(2) = 0.0, x0(3) = -0.2
    nx = 100, ny = 100
 /
```
This input file looks a little more involved than previously.

- We've left the `INPUTPP` section blank as we can simply re-use the
  `methane_charge` file we created when we were generating output for
  xcrysden earlier.
- We've updated `iflag` and `output_format` to select 2D output and gnuplot
  format, and changed the output filename so we'll be able to understand what
  this file was for in future.
- Now since we've selected 2D output, this means we need to define some
  region to plot.
    - The `e1` and `e2` variables are vectors defining a plane
      in units of the lattice constant.
    - `x0` defines the origin of the plane, again in units of the lattice
      constant,
    - `nx` and `ny` are how many points to output in the grid.
    - We've selected an xz plane which contains the carbon and two hydrogens.
      Note, we've shifted the origin slightly so that the carbon will be in
      the centre of the output. You'd need to adjust this to work for the
      shifted structure in the task above.

Now when we run this with `pp.x` we (hopefully) get our `CH4.rho.gpl` data
file as requested. You can take a look at this file, and you'll see it's a
text list of numbers with `x, y, f(x, y)`. So you could use this in any
plotting program or some other software or code for analysis easily.

To plot in gnuplot you can set labels and do a heatmap style plot with the
following gnuplot commands:

```
set title "CH4 Charge Density in y=0 plane"
set xlabel "x Position (Bohr)"
set ylabel "z Position (Bohr)"
set cblabel "Charge Density"
plot "CH4.rho.gpl" with image
```

Note, the x and y coordinates are output in Bohr atomic units as we've set in
the plot here.

Visualizing the Wavefunction
----------------------------

We can also try to visualise the wavefunction we have calculated for
a methane molecule. Take a look at the directory
[`02_wavefunction/01_methane`](02_wavefunction/01_methane). This is contains
an [input file](02_wavefunction/01_methane/CH4.in) for methane exactly as
we've seen before. Note, we've set `disk_io = 'low'` which is the default
value for scf calculations (i.e. we could have omitted this and get the same
output), as we want to keep the wavefunction output file for analysis.

- Run `pw.x` with this input file now, and check the output to make sure it
  all worked as expected. You should have both a `pwscf.save` directory
  containing the charge density file and other output, and a `pwscf.wfc`
  wavefunction file.

Now we want to post-process this output into something we can visualise more
easily. To do that we'll be using the `pp.x` code from the Quantum Espresso
suite. This code can read the output files produced by `pw.x`, extract
whatever quantity you're interested in and generate output compatible with
various visualisation programs. There is a help file similar to the `pw.x`
available in the Quantum Espresso documentation folder called `INPUT_PP.txt`.

### 3D visualization with XCrySDen

There is a `pp.x` input file for methane called
[`CH4_pp.in`](02_wavefunction/01_methane/CH4_pp.in), that contains the
following:

```
 &INPUTPP
    filplot = 'CH4_wfn'
    plot_num = 7
    kpoint(1) = 1
    kband(1) = 1
    kband(4) = 4
 /

 &PLOT
    iflag = 3
    output_format = 5
    fileout = '.xsf'
 /
```

As you can see, the file format is quite similar to that used by `pw.x`. It
contains two sections:

- `INPUTPP` has inputs controlling what data will be extracted from the
  `pw.x` output.
    - If we had chosen non-default directories and names for the output files,
      we would need to set them here.
    - `filplot` as before specifies a filename that will be used to save the
      extracted data.
    - `plot_num` takes an integer option that specifies what quantity will be
      extracted. 7 gives us the wavefunction. There are many options here,
      such as potentials and local density of states. See the help file
      `INPUT_PP.txt` for details, including info on additional variables that
      apply depending on this option.
   - `kpoint(1)` is used to set the k-point index of the wavefunction we want
     to output. There is only the gamma point in our calculation, so we can
     set this to 1. For periodic materials where we have a grid of k-points we
     could also specify `kpoint(2)` and generate output for a range of
     k-points.
   - `kband(1)` and `kband(2)` are the lower and upper limits of the bands
     that we want to output the wavefunction for.
- `PLOT` has inputs controlling how the extracted data will be output for
  visualisation.
   - `iflag` gives the dimensionality of the output. We choose 3 here for a 3D
     plot.
   - `output_format` gives the type of format you want to have output, usually
     this is determined by what application you'll be using to visualise the
     data. It should also be compatible with the value of `iflag`. We've
     selected 5 which will give us output compatible with `xcrysden`.
   - `fileout` usually specifies the name of the plot file, however for
     wavefunction output where several files are generated, the text in
     `fileout` is appended to the name of each file. `xsf` is a default
     `xcrysden` extension, so we can add this here.

We can run `pp.x` in the same way as `pw.x`:

```bash
pp.x < CH4_pp.in &> CH4_pp.out
```

You'll see many files have been generated beginning with `CH4_wfn` followed by
a k-point and band index. Each has a corresponding `.xsf` file. The `.xsf`
files which we'll be using for visualization are in a text format so we can
examine them (and modify them if we want) in a text editor. Take a look and
you'll see the files simply defines the crystal lattice and basis, then have a
section with the datagrid as a mesh of points.

Launch `xcrysden` (recall you'll need to load the xcrysden module first with
`module load xcrysden`), and load the file `CH4_wfn_K001_B001.xsf`. You'll see
a big cube, which has the carbon at each corner and whichever bound hydrogen
falls within the cube nearby. This isn't ideal, but let's carry on for the
moment. There's no sign of a wavefunction yet. To enable this, go to "Tools"
-> "Data Grid". (Make sure you're x2go window and your xcrysdent window are
big enough to see this). You can click "OK" on the menu that appears, and then
you'll see a menu that controls the data grid plot appearance. Try entering
say "0.02" for the "Isovalue" and click "Submit". If you zoom into a carbon
you'll see the wavefunction around it. Unfortunately, since we have the
carbon at the origin, this plot is a little hard to see.

### _Task_

- In a new directory `01_methane_shift`, shift the molecule such that the
carbon atom is at the centre of the box in the `pw.x` input file.
- Run the `pw.x` and `pp.x` calculations as before.
- Visualise the output in `xcrysden` for the wavefunctions associated with
  each of the four bands.
    - Disable the display of the crystal cell.
    - Try plotting the data grid at different isovalues. Each wavefunction
      will have a different range of values that might show interesting
      features. You can also try enabling
    transparency, and increasing the degree of the tricubic cpline to 2 or 3.
    - Under "Modify" try reducing the "Ball Factor" to reduce the visual size
      of the ions.
- When you're happy with how it looks, save an image by first checking the
  "File" -> "Print Setup" menu. (You'll need to leave the data grid menu
  open). It's usually worth enabling anti-aliasing here. Then go to "File" ->
  "Print" and save the output as a png image.

Projected Density of States
---------------------------

Following from the electronic density of states, it can be very useful in
understanding a material to visualize how the density of states can be
decomposed into the various states belonging to each atom in the system. To do
this we can use the `projwfc.x` tool from the Quantum Espresso package. This
is used in a similar way to the `dos.x` tool, although it is not possible to
use tetrahedra to integrate the DOS, so a broadening must be used.

A set of example input files for diamond are given in the directory
[`03_projecteddos/01_diamond`](03_projecteddos/01_diamond). This set of inputs
follows exactly the same progression as previously, except now for the
third step we have an input file for `projwfc.x` rather than for `dos.x`.
In this file we have a single `PROJWFC` section, but we have actually retained
the same inputs we used previously in the DOS calculation. As usual, you can
get full details of the various inputs that are available in the help file
`INPUT_PROJWFC.txt`.

Now run these three calculations as before. Again a short script has been
provided which does this explicitly. Once they have completed you'll see
we again have generated a `pwscf.dos` file as before. You can try to plot this
if you like, and you'll see it's identical to the density of states we
obtained with the equivalent `dos.x` calculation. However, we also have a
number of other files that have been generated, all with names beginning with
`pwscf.pdos_` by default. Try looking at `pwscf.pdos_tot` first. You'll see
this is a three column file, with energy, density of states, and the total of
the various decomposed projected density of states. In principle column 2 and
3 should be the same (and column 2 will reproduce the already calculated
density of states), but in practice it can be difficult to assign conduction
band states accurately so there may be some small disagreement there.

The additional files will the projected density of states for each atom in the
unit cell and each orbital type (s, p, d etc) present. If you look, for
example at the file `pwscf.pdos_atm#1(C)_wfc#1(s)`, you'll see we have three
columns: energy, a column labelled `ldos` and a column labelled `pdos`. For
s-orbitals we only have one value of the magnetic quantum number ml, so there
is only one pdos column, and the ldos and pdos columns are equivalent. If you
look at the file `pwscf.pdos_atm#1(C)_wfc#2(p)` you'll see we now have 5
columns: energy, ldos, and 3 pdos columns. ldos gives the sum of the three
pdos columns, and each pdos column is for a different value of ml (3 for a
p-orbital).

In our case, the corresponding files for each of the C atoms in the cell
should be equivalent as the two atoms in the cell are equivalent. For one
of the atoms, try plotting the projected density of states of the s-orbital
and one of the p-orbitals (e.g. try plotting the 3rd column in each file)
together. Are the states near the top of the valence band more s-like or
p-like?

### _Task_

- Calculate and plot, in whichever way you think best shows the important
  features and differences, the projected density of states for both silicon
  and aluminium.
- How do these compare to each other, and the example diamond calculation?
