Spin Polarization and Magnetic Systems 
======================================

This week we'll cover the topic of spin-polarized systems. For metals, there are
a couple of complications which mean we have to treat them differently from
systems with a non-zero band gap. 

As before, all inpus and scripts you need can be found in
`/opt/MSE404-MM/docs/labs/lab06` and you should make a copy of the folder to
your home directory.

<!-- <div markdown="span" style="margin: 0 auto; text-align: center"> -->
<!-- [Download the input files for this tutorial](./assets/lab08_input.zip){ .md-button .md-button--primary } -->
<!-- </div> -->

------------------------------------------------------------------------------

Spin Polarization
-----------------

Until now we have assumed that the electronic states are the same for up-spin and down-spin electrons: in other words, we assumed that each state can be occupied by two electrons. In a magnetic system, however, the electronic wavefunctions (and the corresponding Kohn-Sham energies) can be different for electrons with different spins. 


The Oxygen Molecule
-------------------

An oxygen molecule has an even number of electrons. So you might expect that a certain number of Kohn-Sham states are filled with two electrons and therefore the total number of up-spin electrons is equal to the total number of down-spin electrons. This is, however, not true. 

When the two oxygen atoms are sufficiently close to each other so that their atomic orbitals start to overlap, new states - called molecular orbitals - are formed. These are called $\sigma_s$, $\sigma_s^*$, $\sigma_z$, $\pi_x$, $\pi_y$ and so on. As a consequence of the symmetry properties that the molecule possesses, some of the molecular orbitals have the same energies. In Quantum Mechanics, such wavefunctions are called degenerate. You can see from the molecular orbital diagram below that the $\pi_x$ and $\pi_y$ orbitals are degenerate and also the $\pi_x^*$ and the $\pi_y^*$ orbitals have the same energy.   

![MO](assets/Mo_diagram.svg){: style="width:250px" align=right}

When we fill these molecular orbitals with electrons, we end up with two electrons that we can distribute in the $\pi_x^*$ and $\pi_y^*$ orbitals. It turns out that the repulsive interaction between electrons favors a state in which the two electrons sit in different molecular orbitals, but have the same spin. This is known as Hund's Rule. Note that a similar situation occurs in an isolated oxygen atom.

The directory `01_O2` contains an input file to calculate the total energy of
the oxygen molecule at the experimentally measured bond length. Here the calculation has been set up
exactly as you've seen in the past (i.e., assuming doubly degenerate band
occupation without smearing or spin polarization):

```python
 &CONTROL
    pseudo_dir = '.'
 /

 &SYSTEM
   ibrav = 1
   A = 10
   nat = 2
   ntyp = 1
   nbnd = 8
   ecutwfc = 60.0
 /

 &ELECTRONS
 /

ATOMIC_SPECIES
 O  15.9999  O.pz-rrkjus.UPF

ATOMIC_POSITIONS angstrom
 O  0.0   0.0  0.0   0 0 0
 O  1.48  0.0  0.0   1 0 0

K_POINTS gamma
```

!!! example "Task 1.1 - Assming Spin Degenerate Insulator"

    Try running the calculation in this directory. Does it converge?
    
    ??? success "Answer"

        While it's possible that the system may randomly meet the convergence
        criteria in the self-consistent cycle, this calculation will most likely
        not converge. If you look at the estimate accuracy at the end of each
        iteration in the output, it will likely vary from step to step, rather
        than steadily decreasing as in a well-behaved calculation.
        
  

To help converge the system, we can use smear the occupancies (as we do in a DFT calculation for a metal). This will allow the system to converge to a ground state.
The relevant input variables are the ones highlighed below:

```python hl_lines="12-14"
 &CONTROL
    pseudo_dir = '.'
 /

 &SYSTEM
   ibrav = 1
   A = 10
   nat = 2
   ntyp = 1
   nbnd = 8
   ecutwfc = 60.0
   occupations = 'smearing' 
   smearing = 'fermi-dirac'
   degauss = 0.1d0
 /

 &ELECTRONS
 /

ATOMIC_SPECIES
 O  15.9999  O.pz-rrkjus.UPF

ATOMIC_POSITIONS angstrom
 O  0.0   0.0  0.0   0 0 0
 O  1.48  0.0  0.0   1 0 0

K_POINTS gamma
```

!!! example "Task 1.2 - Assuming Spin Degenerate Metal"

    Create a copy of the `01_O2` directory called `01_O2_metal`. Modify the
    input file in it to use a metallic occupation scheme with a small smearing
    width and run the calculation (as above). 

    - Does the calculation now converge?

    - Take a look at the file `pwscf.xml` in the calculation directory, and try
      to find the occupations of each band at each k-point. Are these as
      expected?

    ??? success "Answer"
        There are a few states which are completely filled with electrons and then there are two states which are half-filled (i.e. one electron per state). The last state is almost empty.
        ```
                <occupations size="8">
          9.999997613058770E-001  9.999680732750561E-001  9.800308333633008E-001
          9.433101748955524E-001  9.433101708179502E-001
          5.459533880551336E-001  5.459533513549418E-001  4.147424691420094E-002
                </occupations>
        ```

While treating this system as a metal may help converging the calculation, it
may not necessarily reach the true ground state (i.e. the state with the lowest total energy) since we have imposed that there is an equal number of up-spin and down-spin electrons. To allow different numbers of spin-up and spin-down electrons, we can perform a spin-polarized calculation by adding the
`nspin` and `tot_magnetization` variables to the input file (highlighted below):

```python hl_lines="12-13"
 &CONTROL
    pseudo_dir = '.'
 /

 &SYSTEM
   ibrav = 1
   A = 10
   nat = 2
   ntyp = 1
   nbnd = 8
   ecutwfc = 60.0
   nspin = 2 #(1)!
   tot_magnetization = 2.0 #(2)!
 /

 &ELECTRONS
 /

ATOMIC_SPECIES
 O  15.9999  O.pz-rrkjus.UPF

ATOMIC_POSITIONS angstrom
 O  0.0   0.0  0.0   0 0 0
 O  1.48  0.0  0.0   1 0 0

K_POINTS gamma
```

1.    `nspin`: this is 1 by default so no spin polarization is taken into 
      account. To perform a spin polarized calculation it should be set to 2.
2.    `tot_magnetization`: this is difference between the number of spin-up and
      spin-down electrons in the cell.

!!! example "Task 1.3 - Assuming Spin Polarized Metal"

    Create another copy of `01_O2` called `01_O2_spin`. Then, try to:

    1. Only turn on spin polarization (`nspin=2`). Does the calculation run?

        ??? success "Answer"
            The calculation will not run because it needs to know how to set the
            number of electrons in each spin channels.
            ```
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Error in routine iosys (1):
                fixed occupations and lsda need tot_magnetization
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ```

    2. Setting the total magnetization to 0.0, which would be the case if we 
       don't have any net magnetization in the molecule, as both spins point in
       opposite directions. What is the final energy?

        ??? success "Final energy"
            The calculation converges to an energy of -63.25520699 Ry.

    3. Setting the total magnetization to 2.0, which corresponds
       to both spins pointing in the same direction. Is the energy lower? How do
       the orbital energies vary?

        ??? success "Answer"
            The calculation converges to an energy of -63.29338911 Ry. The
            energy becomes lower with this configuration and the orbital
            energies become different between spin channels. Note that the
            eigenvalues and occupations are written in the `pwscf.xml` file and
            the spin up values are written first, followed by the spin down 
            values.
            ```
                  <ks_energies>
                    <k_point weight="1.00000000000000">0.000000000000000E+000  0.000000000000000E+000  0.000000000000000E+000</k_point>
                    <npw>26462</npw>
                    <eigenvalues size="16">
             -1.050830548259423E+000 -8.119807870761847E-001 -4.770004419270685E-001
             -4.341811684674406E-001 -4.341811560749695E-001
             -3.082114385895908E-001 -3.082114249408998E-001 -1.280528718472738E-001
             -9.970063621108551E-001 -7.461771608702750E-001
             -4.381689137303877E-001 -3.691144532629034E-001 -3.691144467746303E-001
             -2.318554649322491E-001 -2.318554571437267E-001
             -8.403816145565962E-002
                    </eigenvalues>
                    <occupations size="16">
              1.000000000000000E+000  1.000000000000000E+000  1.000000000000000E+000
              1.000000000000000E+000  1.000000000000000E+000
              1.000000000000000E+000  1.000000000000000E+000  0.000000000000000E+000
              1.000000000000000E+000  1.000000000000000E+000
              1.000000000000000E+000  1.000000000000000E+000  1.000000000000000E+000
              0.000000000000000E+000  0.000000000000000E+000
              0.000000000000000E+000
                    </occupations>
                  </ks_energies>
            ```

Finally, comparing the energy of the spin polarized calculation with the spin
degnerate metal calculation, we can see that the spin polarized calculation
gives a lower energy.

!!! pied-piper "Fun fact"

    O2 in its singlet state can be dangerous (see e.g. [:link:`this
    paper`](https://www.sciencedirect.com/science/article/pii/S1383574211001189)),
    so treating the spin correctly is important!


Iron :material-hammer:
----------------------

Now that you've seen how including spin polarization can allow us a correctly
describe the ground state of a molecular system, the next step is to use it to
describe a magnetic crystal system.

In a magnetic crystal there is a net spin polarization in the unit cell. This
means that we'll probably have an odd number of electrons, and the energy of the
system when we include a net spin polarization is lower than the energy when we
don't.

One of the most common magnetic crystal is iron, so we'll examine this. The
directory `02_Fe` contains an input file for iron. Note this is a BCC structure
(as set by `ibrav = 3` in the input file), whereas most of the crystals
structures you have examined previously were FCC. The calculation has been set
up in the usual way for a metallic system.

!!! example "Task 2.1 - Fixed Magnetization"

    1. Run this calculation and check everything worked as expected. What is the
       final energy?

        ??? success "Answer"
            The final energy should be -55.52528610 Ry.

    2. Now make a copy of the calculation directory and in this, modify the
       calculation to turn on spin polarization. Try running the calculation
       with `tot_magnetization = 0.0` first, and compare your total energy to
       that obtained using doubly degenerate bands. 

        !!! Warning "Warning" 
            While in the case of the O2 above, we were able to get our
            calculations to at least converge by using a metallic occupation
            instead of using spin polarization, in the case of iron, it will
            still be a metal when you use spin polarization, so you should not
            remove the input variables associated with this. 

        ??? success "Answer"
            The total energy becomes -55.52528589 Ry. Almost identical to the
            one obtained with the doubly degenerate bands. This is because these
            two calculations are essentially identical.

    3. Now try setting the total magnetization to 1.0 and see how total energy
       changes: Which is the more energetically favourable configuration?

        ??? success "Answer"
            The total energy becomes -55.53839616 Ry. Lower than the spin
            degenerate case.

    4. Try setting the total magnetization to 2.0. How does the final energy
       compare to the previous value?

        ??? success "Answer"
            The total energy becomes -55.56226730 Ry. Lower than all previous
            cases.

From this we could test many guesses for the total magnetization, and find
the value which gives the lowest overall total energy. Note that here one can
set the total magnetization to be a fractional number which is not physical for
molecules. This is because we have now a periodic metal systems where
itinerant electrons can also be a media to host spin polarization. 

However, finding the ground state total magnetization value can be a trdious job
and one can instead pass an option that tells quantum espresso to automatically 
find the best value. This is done by setting the `starting_magnetization` input
variable.

!!! example "Task 2.2 - Relaxed magnetization"

    1. Make another copy of the `02_Fe` directory, and this time set `nspin =
       2`, and `starting_magnetization = 1.0` (do not include the
       `tot_magnetization` variable as this fixes a value). Run the calculation
       and see what the final total magnetization per cell is. See if you can
       find a measured value for iron to compare to.

        ??? success "Answer"
            The total magnetization becomes larger than 2.0.
            ```
            total magnetization       =     2.21 Bohr mag/cell
            ```
            This is becuase we are allowing the spin to fully relax in the
            system.
    
    2. See if you can use what we covered in previous labs to calculate and make
       a plot of the electronic band structure of BCC Fe.

        - Plot the spin-up and spin-down bands in different colours.
        - Indicate the Fermi energy on your plot in some sensible way.
        - As the Brillouin zone is different to the ones you have calculated so
          far you'll need to select a few sensible high-symmetry points yourself
          to plot with :slight_smile:.

        ??? success "Answer"
            You can find the relevant input file in the directory
            `02_Fe/extra_bandstructure`. The band structure should look similar
            to the following: 
            <figure markdown="span"> 
            ![Diamond primitive cell](assets/Iron_bands.png){ width="500" } 
            </figure>
            You can find `README.md` in the directory for more information on
            how to reproduce this plot.

------------------------------------------------------------------------------

Summary
-------

In this lab you have seen:

<!-- - How to treat a metallic system. -->
- How to do a DFT calculation including spin polarization.
- How some systems need to be done with spin polarization to converge to the
  correct ground state.
- How to use spin polarized calculations to find the correct magnetization of a
  magnetic system by letting the code find the total magnetization that produces
  the lowest overall total energy.
