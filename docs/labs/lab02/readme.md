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
        ibrav is now set to 0. This means 'free cell' meaning the user needs to specify the cell parameters manually. The parameter defining the cell vector length is not present. There is also a section called `CELL_PARAMETERS`.

    - Will this input file do the exact same thing as the one in `01_methane`?

    ??? success "Answer"
        Yes! Instead of specifying the cell vector length we have just specified the length of each cell vector in the `CELL_PARAMETERS` section. 


!!! warning "Warning - Periodic Boundary Conditions and Molecules"
    Quantum Espresso is a periodic DFT code due to it using a plane wave basis. This is something you will learn about later in the theoretical part of this course. We therefore need to be smart in order to model isolated molecules, since by definition these are not periodic. One way of doing this is to put the molecule in the center of a large box. This minimises any interaction with its periodic neighbour.


## Running and examining the calculation

The Quantum Espresso package has been compiled as a module on the server as discussed in [Lab 1](../lab01/readme.md). Modules are often used on HPC systems to make different versions of packages available to users.
In order to be able to run any Quantum Espresso calculation, it must first be loaded to your environment. This can be done by issuing the command

```bash
module load quantum-espresso
```

This will load Quantum Espresso and any module dependencies.

!!! example "Task 3 - Running a calculation"
    To run the first calculation of the day, make sure you have loaded Quantum Espresso to your environment as discussed above.

    - Navigate to the `01_methane` directory.

    - Issue the command 
    ```bash
    pw.x < CH4.in > CH4.out
    ```
    
After the calculation has finished take a look at the files created in your directory. You should have a file named `pwscf.xml` and a new directory named `pwscf.save`.

- `pwscf.xml` has the results of the pw.x calculation in machine readable format (not so readable for humans!)
- `pwscf.save` is a directory which contains: 
	- A copy of `pwscf.xml`.
	- A copy of the pseudopotential files used in the calculation.
	- A file with the charge density stored. This is mostly used for post-processing to obtain observables.
	- Wavefunction files which are stored in binary (and thus are not readable). These can be used as inputs to other calculations.

Now that we have run the calculation for methane, we should examine the output file `CH4.out`, which is where we instructed Quantum Espresso to pipe the output of the pw.x calculation.
Using the command
```bash
more CH4.out
```
we can look at the output file. Output files are generally structured as such:
- Begining - Important information about the system including parameters used in the calculation.
- Middle - Self-consistent cycle starting from randomised atomic wavefunctions.
- End - Final results like total energy, band energies etc.

!!! example "Task 4 - Examining an output file"
    Using the `more` command specified above

    - How many electrons were in your calculation?

    ??? success "Answer"
        8.00. This is found at the top of your output file.
        ```bash
        number of electrons       =         8.00
        ```

    - How many scf cycles (iterations) did your calculation go through?

    ??? success "Answer"
        9 scf cycles. This is found on the line:
        ```bash
        convergence has been achieved in   9 iterations
        ```

    - What is the total energy of the Methane molecule?

    ??? success "Answer"
        $E_{\text{Tot}} = -15.49834173 \, \text{Ry}$. This is found on the line:
        ```bash
        !    total energy              =     -15.49834173 Ry
        ```
        Notice the converged total enegry will always have a `!` at the beginning of the line.

    - What accuracy is your calculation converged to?

    ??? success "Answer"
        0.00000066 Ry. This is found on the line:
        ```bash
        estimated scf accuracy    <       0.00000066 Ry
        ```
        We did not specify this in the input file. The default value of below 1E-6 was therefore used.

    - How many band energies were calculated?

    ??? success "Answer"
        4 band energies were calculated. This is found in the lines:
        ```bash
        End of self-consistent calculation

        k = 0.0000 0.0000 0.0000 ( 14712 PWs)   bands (ev):

        -17.3307  -9.3182  -9.3176  -9.3173
        ```

??? note "Electrons and Energy Eigenvalues"
    Note that in this calculation we had 8 electrons but only 4 energy eigenvalues were calculated. This is because we have treated the 8 electrons as 4 doubly occupoed states, and therefore only 4 energy eigenvalues are outputted.

!!! example "Task 5 - Alternative Input File"
    Navigate to the directory `01a_methane`. 

    - Run the same calculation as in Task 4 and confirm that you get the same results.

## Visualising Structures - VESTA

Interactive visualisation software are highly important in computational physics. Not only are they a way of checking the structure defined in your input file, but they are also very useful when checking output structures of relaxation calculation. You will learn more about this in [Lab 5](../lab05/readme.md). The visualisation software we are going to use through this course is called `VESTA`.

Vesta, like Quantum Espresso, has been loaded into a module. To use it you will need to issue the command:

```bash
module load vesta
```

You have now loaded VESTA to your environment. By default, VESTA cannot read Quantum Espresso input files. Therefore, we will need to convert to a format that VESTA can read.

***I don't yet know how we are going to convert the QE to VESTA file.***

After converting your file, you will be able to visualise it with the command:

```bash
vesta filename.whatgoeshere???
```

During this lab we will be working with different molecules. It will be a good exercise to visualise them as we go along.

## Methane, ethane and ethene

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




