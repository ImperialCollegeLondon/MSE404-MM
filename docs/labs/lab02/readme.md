Quantum Espresso Input and Output for Molecules
===============================================

This week you will run some DFT calculating using Quantum Espresso for molecular systems. We are going to be focusing on understanding and constructing input files that you will need for the rest of the course.

<div markdown="span" style="margin: 0 auto; text-align: center">
[Download the input files for this tutorial](./assets/lab02_input.zip){ .md-button .md-button--primary }
</div>

Before starting, if you can't remember how to do something from the command line, you can always refer back to [lab 1](../lab01/readme.md).

------------------------------------------------------------------------------

## Quantum Espresso

[Quantum Espresso](http://www.quantum-espresso.org) is a freely available package of open-source codes for electronic-structure calculations and materials modelling at the nanoscale. It is based on density-functional theory, plane waves, and pseudopotentials, which you will be learning about in lectures.

Quantum Espresso is used via the command line. There is no graphical interface by default, which is typical of most electronic structure codes. Throughout this course we will be interacting with Quantum Espresso through the command line via ssh (through PuTTY that you used in [lab 1](../lab01/readme.md)).


!!! example "Task 1 - Copy Input Files"

    In lab 1 you should have created a directory named `MSE404` in your home directory.

    - Check this by issuing the command `cd ~` followed by `ls`.
    - Copy the input files from `opt/Courses/MSE404/lab02` to your `MSE404` folder. Remember you need to pass an additional flag to `cp` to copy a directory. If you are struggling with this, revisit [lab 1](../lab01/readme.md).
    - Copy the directory containing the pseudopotentials that you will be using during this course to your `MSE404` directory. These are stored in `/opt/Courses/MSE404/pseudo`

You should now have a directory `lab02` and `pseudo` within your `MSE404` directory. This contains a set of basic input files for a variety of systems and the pseudopotentials for the input files.

## Input Files

Before running a calculation we need to write input files. These give instructions to Quantum Espresso to tell it what we want to calculate, and what parameters to use to do the calculation. The first example we will be looking at is in the [01_methane](01_methane) directory. This is an input file for the `pw.x` module of Quantum Espresso which calculates the total energy of your system.

Let's take a look at our first input file [CH4.in](01_methane/CH$.in).

```bash
&CONTROL #(1)!
    pseudo_dir = '.' #(2)!
/

&SYSTEM
    ibrav = 1 #(3)!
    A = 15.0 #(4)!
    nat = 5 #(5)!
    ntyp = 2 #(6)!
    ecutwfc = 18.0 #(7)!
/

&ELECTRONS
/

ATOMIC_SPECIES
 C  12.011  C.pz-vbc.UPF #(8)!
 H   1.008  H.pz-vbc.UPF

ATOMIC_POSITIONS angstrom #(9)!
 C  0.0       0.0       0.0
 H  0.0       0.0       1.089
 H  1.026719  0.000000 -0.363000
 H -0.513360 -0.889165 -0.363000
 H -0.513360  0.889165 -0.363000

K_POINTS gamma #(10)!
```

1. Quantum Espresso input files are ordered with 'tags'. These 'tags' start with a `&` and end with a `/`. They are blocks of input parameters.
2. Directory containing your pseudopotentials defined later in the input file. The directory `.` means the current directory. 
3. Bravais lattice type e.g. simple cubic, face centered cubic etc. These are documented on the [Quantum Espresso input description page](https://www.quantum-espresso.org/documentation/input-data-description/). You will get familiar with this throughout the course. ibrav = 1 is a simple cubic lattice.
4. Crystalographic constant i.e. cell vector length. Simple cubic with A=15 means a 15x15x15 Å box. 
5. Number of atoms.
6. Number of species.
7. Energy cutoff for wavefunction expansion. You will learn more about this in your lectures and [lab03](../lab03/readme.md)
8. Atomic species, atomic mass and the name of the pseudopotential file.
9. The atomic positions of your atoms. The `angstrom` after `ATOMIC_POSITIONS` means these are cartesian coordinated in units of Å.
10. K-Points for the calculation. We have chosen to do this calculation at the Gamma point.

!!! example "Task 2 - Alternative Input File Style"
    Take a look at the input file in the `01a_methane` directory.

    - How is this different from the input file discussed above?

    ??? success "Answer"
        The parameter defining the cell vector length is not present. There is also a section called `CELL_PARAMETERS`.

    - Will this input file do the exact same thing as the one in `01_methane`?

    ??? success "Answer"
        Yes! Instead of specifying the cell vector length we have just specified the length of each cell vector in the `CELL_PARAMETERS` section. 


!!! warning "Warning - Periodic Boundary Conditions and Molecules"
    Quantum Espresso is a periodic DFT code due to it using a plane wave basis. This is something you will learn about later in the theoretical part of this course. We therefore need to be smart in order to model isolated molecules, since by definition these are not periodic. One way of doing this is to put the molecule in the center of a large box. This minimises any interaction with its periodic neighbour.


## Running the calculation


The Quantum Espresso package has been compiled as a module on the mt-student
server. As discussed in the previous lab, modules are often used on HPC
systems to make different versions of various packages as compiled with
different compilers available to users. To add Quantum Espresso to your
environment, along with its dependencies type the following in a terminal:

```bash
module load gcc mkl espresso
```

Now to run the code: make sure you are in the `01_methane` directory
then do:

```bash
pw.x < CH4.in
```

Here we redirected our input file to the stdin of `pw.x`, as described in the
[IO Redirection section of lab 1](../lab01/readme.md#io-redirection). (We could
also have used the `-i` flag to specify an input file as `pw.x -i CH4.in`.)
You'll see a lot of output has been generated in your terminal, but a simple
calculation like this will finish quite quickly. It's better to explicitly save
the output (and any error information) to a file. To do this, we can instead run
the calculation a redirect the output to a file with

```bash
pw.x < CH4.in &> CH4.out
```

The output files
----------------

Take a look in the directory with `ls`. You'll see some additional files have
been generated. We have `pwscf.xml`, which has all the details of the system
and results of the calculation in an xml format. If you skip to the end of
this file (if you have opened it with `less`, then pressing `G` will go
straight to the end, while `g` will go back to the start) you will see the
eigenvalues and other details. _Note, that while the band energies listed in the
output file are in eV, those in the xml file are in Hartree units_. And we have
a folder called `pwscf.save`. This has some other files such as the charge
density, another copy of the xml file we saw above, a copy of the pseudoptential
that was used, and files with the calculated wavefunction stored in a binary
format that can be used as input to other calculations.

Now let's take a quick look through the output that we generated. (e.g.
`less CH4.out`).

- First we have details of the code, version and date of the calculation.
- There's some info about the calculation details, including some numbers
  where we had just accepted the defaults.
- Then we go into our self-consistent calculation, starting from randomized
  atomic wavefunctions. At the end of each step it outputs the total energy.
- After seven steps we have converged to the default level, i.e. 
  `estimated scf accuracy` is below 1.0E-6.
- The eigenvalues are output in eV. Note although we have
  8 electrons, we have treated them as 4 doubly occupied states, so only four
  numbers are output.
- The total energy is listed, as well as its decomposition into several terms.
- And finally some timing information. If calculations are taking a long time
  this can be helpful in seeing where it is spent.

Note there are a number of equivalent ways of defining the unit cell vectors
within Quantum Espresso. Instead of relying on choosing the correct `ibrav`
value, we can set `ibrav = 0` and give the unit cell lattice vectors
explicitly in a `CELL_PARAMETERS` section. This is shown in
[`01a_methane/CH4.in`](01a_methane/CH4.in).

### _Task_

- Try running `pw.x` on the second methane input file and confirm you get the
  same total energy as before.


Methane, ethane and ethene
---------------------------

Now we understand the basics of the Quantum Espresso input file, let's try
some other molecules, in this case ethane and ethene. The only things that we
need to change are the number of atoms and the atomic positions. The input 
files are in [C2H6.in](02_ethane/C2H6.in) and [C2H4.in](03_ethene/C2H4.in). If
you want to see what we've changed in the [ethane input file](02_ethane/C2H6.in)
relative to the [methane input file](01_methane/CH4.in) a useful tool is
[`diff`](../extras/misc/linuxcommands/readme.md#diff). If you're in the `lab02`
directory you can type `diff 01_methane/CH4.in 02_ethane/C2H6.in`. You'll be
able to see that we've only changed a few lines in our input file and everything
else is the same.

### _Task_

- Run `pw.x` for ethane and ethene and save the outputs.
- How do the eigenvalues compare between molecules?

Note that, unlike for methane, we specified the atomic units in bohr for both 
ethane and ethene. One common mistake is using the wrong units. What happens if
you get the units wrong?

### _Task_

- Change the units of the atomic positions for ethene and run
  `pw.x` again. What happens?


C<sub>20</sub> isomers
----------------------

While the total energy of a molecule isn't that useful by itself, the _relative_
energies between, say, different isomers of a given molecule are much more useful.
In general (ignoring e.g. temperature effects), a lower energy indicates a more
stable isomer. As an example, let's consider three different isomers of
C<sub>20</sub> - a bowl, ring and cage structure
(see <https://pubs.acs.org/doi/abs/10.1021/acs.jpca.5b10266>).

### _Task_

- Run the inputs for the different isomers, found in
  [`C20_bowl.in`](04_c20_bowl/C20_bowl.in), [`C20_ring.in`](05_c20_ring/C20_ring.in)
  and [`C20_cage.in`](06_c20_cage/C20_cage.in). Which one has the lowest total energy?

- The three isomers are close in energy, so how can we tell if the difference is
  significant? Let's try looking at an amorphous structure.

- Run the inputs for the amorphous structure, found in
  [`C20_amorphous.in`](07_c20_amorphous/C20_amorphous.in). How does the energy
  compare to the previous isomers? Is it what you expected?

- Finally, let's see what happens if we try something less realistic. In 
  [`C20_smile.in`](08_c20_smile/C20_smile.in) there is an input file for a structure
  which resembles a smiley face. Run this calculation and compare the energies.

It turns out that the relative energies of the ring, cage and bowl structure are very
sensitive to the details of the approach used, as seen in
[this paper](https://pubs.acs.org/doi/abs/10.1021/acs.jpca.5b10266).
For these systems, it is therefore hard to say with certainty which is the most
stable isomer using DFT. However, by comparing with less realistic structures we can
still see the difference in energy between structures which should be stable
and those which should not.

Finally, note that if you were to use a structure which more closely resembles a smiley
face, Quantum Espresso is unable to reach convergence within 100 iterations, just like
when you use the wrong units. So if you encounter convergence problems, the first thing
to check should be whether or not your input structure is sensible.

Visualizing Structures
----------------------

There are many different tools that can be used to visualize atomic
structures. `xcrysden` is installed as a module on the server you're using for
this course, and conveniently can read Quantum Espresso input files. Try
loading the module with `module load xcrysden`, running the command `xcrysden`
and opening the input files for the various structures we've looked at in this
lab. You can do this by selecting `Open PWscf` on the file menu. Note that
since we are viewing a molecule, xcrysden will open a menu asking whether to
reduce the dimensionality. Select `reduce dimension to 0D`. 
There are many options to control how the structure looks, and you can
grab and rotate the structure with your mouse.

### _Task_

- See if you can figure out how to save an image of each C<sub>20</sub> isomer
  as a png file.

We can also use `xcrysden` to visualize other quantities, such as the charge
density. If you've got time left in this lab, check out the additional material
on [`visualization`](../extras/labs/visualising_output/readme.md), which shows
you how to visualize the charge density of methane.

-------------------------------------------------------------------------------

Summary
-------

- In this lab we've looked at how to create an input file for some different
  molecules:
  - Methane, ethane and ethene.
  - Different isomers of C<sub>20</sub>.
- We've also looked at how to visualize the structures.
- As an optional extra, try setting up inputs for some other molecules of your
  choice.


-------------------------------------------------------------------------------

Extra - Other molecules
-----------------------

In this course we give you the coordinates of the materials that you need to 
simulate. But what if you want to try something different? For molecules,
[Avogadro](https://avogadro.cc/) is a useful for progam for visualizing
and generating structures for a range of molecules. The online documentation
will show you how to build different molecules. Once you are happy, you can 
save the coordinates as a `.xyz` file, which will write the atomic positions
in Angstrom. You can then insert these coordinates into a Quantum Espresso
input file.

### _Optional Task_

- Try building a (small) molecule of your choice in Avogadro. For now, just use
  carbon and hydrogen - we'll look at what we need to do to use different elements
  next week. Once you've built a molecule, use it to make a Quantum Espresso input
  file and run the calculation.




