Spin Polarization and Magnetic Systems 
======================================

This week we'll cover the topic of spin-polarized systems. For metals, there are
a couple of complications which mean we have to treat them differently from
systems with a non-zero band gap.

<div markdown="span" style="margin: 0 auto; text-align: center">
[Download the input files for this tutorial](./assets/lab08_input.zip){ .md-button .md-button--primary }
</div>

------------------------------------------------------------------------------

<!-- ## Metals -->
<!--  -->
<!-- Metals have a Fermi surface that can be quite complex in k-space. This means -->
<!-- that in contrast to an insulator or semiconductor where every k-point has the -->
<!-- same number of occupied states, **in a metal the number of occupied states can -->
<!-- vary from k-point to k-point**. This makes them more difficult to converge than -->
<!-- other systems.  -->
<!--  -->
<!-- In short, there are generally two things you need to do: -->
<!--  -->
<!-- 1.  Use a denser k-point grid than you would need for a semiconductor or -->
<!--     insulator. This is to help sampling the rapid change in the Fermi surface at -->
<!--     different k-points. -->
<!--  -->
<!-- 2.  Use some smearing scheme. This is in relation to the smearing used in the -->
<!--     calculation of the [:link: density of -->
<!--     states](../lab04/readme.md#density-of-states). The difference is that here -->
<!--     the occupation is also smeared (i.e., can no longer be intergers of 0 and -->
<!--     1.) Visually, the smeared DOS would look like the following: -->
<!--  -->
<!--     where the occupation function (Fermi-Dirac function) is plotted in red. The -->
<!--     Fermi energy is obtaeind by solving the following equaion: -->
<!--     $$  -->
<!--     \int_{-\infty}^{\varepsilon_F} \mathrm{DOS}(\varepsilon) f_T(\varepsilon) -->
<!--     d\varepsilon = N_e  -->
<!--     $$ -->
<!--     where $N_e$ is the number of electrons in the system and $f$ represents the -->
<!--     Fermi-Dirac distribution function at temperature $T$. As we already know,  -->
<!--     the Fermi-Dirac function at 0K is a step function which would spoil the -->
<!--     convergence of metals. Here, we simply raise the temperature to a small -->
<!--     number (`degauss`) so that the Fermi-Dirac function is smeared out and the -->
<!--     Convergence can be more easily achieved. It is worth noting that other -->
<!--     smearing methods such as gaussian smearing can also be used. -->
<!--  -->
<!--     Adding a smearing helps significantly in achieving a smooth -->
<!--     SCF convergence, as otherwise a small change in a state energy from once -->
<!--     cycle to the next could lead to a very large change in its occupation and to -->
<!--     the total energy in turn (this is called 'ill-conditioning').  -->
<!--     We set the smearing scheme and width with the `occupations` and `degauss`  -->
<!--     variables in the input file. -->
<!--  -->
<!-- Example: Aluminium -->
<!-- ------------------ -->
<!--  -->
<!-- Aluminium forms in a standard fcc structure with one atom per cell, which we -->
<!-- know how to deal with at this point. The thing about Aluminium that makes it -->
<!-- more complicated within DFT is that it is a metal. -->
<!--  -->
<!-- Here is an example input file for a calculation of Aluminium: -->
<!--  -->
<!-- ```python -->
<!--  &CONTROL -->
<!--     pseudo_dir = '.' -->
<!--  / -->
<!--  -->
<!--  &SYSTEM -->
<!--     ibrav =  2 -->
<!--     A = 2.863 -->
<!--     nat =  1 -->
<!--     ntyp = 1 -->
<!--     ecutwfc = 18.0 -->
<!--     occupations = 'smearing' #(1)! -->
<!--     smearing = 'fermi-dirac' #(2)! -->
<!--     degauss = 0.1d0 #(3)! -->
<!--  / -->
<!--  -->
<!--  &ELECTRONS -->
<!--  / -->
<!--  -->
<!-- ATOMIC_SPECIES -->
<!--  Al  26.982  Al.pz-vbc.UPF -->
<!--  -->
<!-- ATOMIC_POSITIONS crystal -->
<!--  Al 0.00 0.00 0.00 -->
<!--  -->
<!-- K_POINTS automatic -->
<!--   8 8 8 1 1 1 -->
<!-- ``` -->
<!--  -->
<!-- 1.    The `occupations` variable is set to `smearing` to tell Quantum Espresso -->
<!--       to use a smearing scheme [:link:input  -->
<!--       description](https://www.quantum-espresso.org/Doc/INPUT_PW.html#idm362). -->
<!-- 2.    The `smearing` variable is set to `fermi-dirac` to tell Quantum Espresso -->
<!--       to use a Fermi-Dirac smearing scheme. [:link:input -->
<!--       description](https://www.quantum-espresso.org/Doc/INPUT_PW.html#idm404).  -->
<!-- 3.    The `degauss` variable is set to 0.1d0 to set the width of the smearing. -->
<!--       see [:link:input -->
<!--       description](https://www.quantum-espresso.org/Doc/INPUT_PW.html#idm401). -->
<!--  -->
<!--  -->
<!-- !!! example "Task 1 - Smearing" -->
<!--  -->
<!--     First, run the `pw.x` calculation with the supplied input file in -->
<!--     [:link:01_aluminium/Al.in](01_aluminium/Al.in). -->
<!--      -->
<!--     Then, look in the `pwscf.xml` file and find the various `ks_energies` -->
<!--     entries towards the end. These give the various k-points used in the -->
<!--     calculation and the energies and occupations of each state for this k-point. -->
<!--     Note, for a metal the default number of bands is at least four more than are -->
<!--     needed for the number of electrons per cell. The pseudopotential we have -->
<!--     used has 3 valence electrons, which could be represented with two -->
<!--     potentially doubly occupied bands, so we have four more bands in the -->
<!--     calculation for a total of 6. -->
<!--  -->
<!--     ??? success "Example"  -->
<!--         ``` -->
<!--               <ks_energies> -->
<!--                 <k_point weight="7.812500000000000E-003">-6.250000000000000E-002  6.250000000000000E-002  6.250000000000000E-002</k_point> -->
<!--                 <npw>59</npw> -->
<!--                 <eigenvalues size="6"> -->
<!--           1.315972343567215E-001  1.505697520824042E+000  1.607697079464305E+000 -->
<!--           1.607697714947740E+000  1.834366371282428E+000 -->
<!--           1.952726961146777E+000 -->
<!--                 </eigenvalues> -->
<!--                 <occupations size="6"> -->
<!--           9.999990177787399E-001  1.181697427742303E-006  1.536561074875367E-007 -->
<!--           1.536541545820267E-007  1.650917762173208E-009 -->
<!--           1.547598926179030E-010 -->
<!--                 </occupations> -->
<!--               </ks_energies> -->
<!--          ...  -->
<!--         ``` -->
<!--      -->
<!--     Now, try removing the `occupations` and `degauss` variables from the input -->
<!--     file and see what happens when you try to run the calculation. -->
<!--  -->
<!--     ??? success "Example"  -->
<!--         ``` -->
<!--         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -->
<!--              Error in routine electrons (1): -->
<!--              charge is wrong: smearing is needed -->
<!--         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -->
<!--         ``` -->
<!--  -->
<!-- ------------------------------------------------------------------------------ -->

Spin Polarization
-----------------

Up untill now we have been assuming that we always had some set of bands which
could each fit two electrons. Essentially we have been ignoring the electron
spin. If you want to examine, for example, a magnetic system then the spin of
the electrons is important. It can also be important in modelling atomic or
molecular systems. We'll cover different examples of this in this lab.


The Oxygen Molecule
-------------------

If a system is not necessarily magnetic we might imagine that representing 
it with some set of fully occupied, doubly degenerate bands will work. However,
in some cases including spin polarization can lead to important differences. One
example of this is the O2 molecule.

![MO](assets/Mo_diagram.svg){: style="width:250px" align=right}

In this case, we have a system with two interacting oxygen atoms. Each oxygen
has 8 electrons in total, with the configuration `1s2 2s2 2p4` (the 1s orbital
will be contained within the pseudopotential for the DFT calculations done
here, so you will have 6 electrons from each oxygen atom). 

For a single oxygen, from Hund's rule the three p orbitals should be filled
singly before being filled in pairs, so that one of the p-orbitals will have two
electrons, and the other two should have one each. However, if we assume doubly
occupied orbitals, we'll have the two p-orbitals with two electrons and one that
is empty. This means a calculation where we assume a set of doubly occupied
bands will have trouble converging to the ground state of the system. For the
molecule the situation is similar, but the s and p orbitals from each atom
combine to form bonding and anti-bonding $\sigma$ and $\pi$ orbitals.

The directory `02_O2` contains an input file to calculate the total energy of
the system at the measured bond length. Here the calculation has been set up
exactly as you've seen in the past (i.e., assuming doubly degenerate band
occupation without smearing or spin polarization:

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

!!! example "Task 2.1 - Assming Spin Degenerate Insulator"

    Try running the calculation in this directory. Does it converge?
    
    ??? success "Answer"
        While it's possible that the system may randomly meet the convergence
        criteria in the self-consistent cycle, this calculation will most likely
        not converge. If you look at the estimate accuracy at the end of each
        iteration in the output, it will likely vary from step to step, rather than
        steadily decreasing as in a well-behaved calculation.
        
        The situation we have is similar to a metal: we have two bands and the ground
        state of the system should be when there is one electron in each of them.

To get around this, we can use a metallic occupation scheme with a small
smearing width. This will allow the system to converge to the correct ground
state. The relevant input variables are the ones highlighed below:

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

!!! example "Task 2.2 - Assuming Spin Degenerate Metal"

    Create a copy of the `01_O2` directory called `01_O2_metal`. Modify the
    input file in it to use a metallic occupation scheme with a small smearing
    width and run the calculation (as above). 

    - Does the calculation now converge?

    ??? success "Answer"
        Yes, the calculation should now converge.

    - Take a look at the file `pwscf.xml` in the calculation directory, and
      try to find the occupations of each band at each k-point. Are these as
      expected?

    ??? success "Answer"
        The occupations should be fractional for the highest occupied valence
        band which is not physical for a molecule.
        ```
                <occupations size="8">
          9.999997613058770E-001  9.999680732750561E-001  9.800308333633008E-001
          9.433101748955524E-001  9.433101708179502E-001
          5.459533880551336E-001  5.459533513549418E-001  4.147424691420094E-002
                </occupations>
        ```

While treating this system as a metal may help converging the calculation, it
may not necessarily reach the ground state since the spin-degress of freedom is
constrained. Instead, we can do a spin polarized calculation by adding `nspin`
and `tot_magnetization` variables to the input file (highlighted below):

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
      spin-down electrons in the cell. If we want a single spin up electron
      we can set this to `1.0`.

!!! example "Task 2.3 - Assuming Spin Polarized Metal"

    Create another copy of `02_O2` called `02_O2_spin`. Then, try to:

    1. Only turn on spin polarization. Does the calculation run?

        ??? success "Answer"
            The calculation will not run.
            ```
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Error in routine iosys (1):
                fixed occupations and lsda need tot_magnetization
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ```

    2. Setting the total magnetization to 0, which would be the case if we 
       don't have any net magnetization in the molecule, as both spins point in
       opposite directions. 

        ??? success "Answer"
            The calculation converges to an energy of -63.25520699 Ry.

    3. Setting the total magnetization to 2.0, which corresponds
       to both spins pointing in the same direction. Is the energy lower? How do
       the orbital energies vary?

        ??? success "Answer"
            The calculation converges to an energy of -63.29338911 Ry. The
            energy becomes lower with this configuration and the orbital
            energies become different between spin channels.
            ```
                  <ks_energies>
                    <k_point weight="1.00000000000000">0.000000000000000E+000  0.000000000000000E+000  0.000000000000000E+000</k_point>
                    <npw>26462</npw>
                    <eigenvalues size="16">
             -1.025922874232402E+000 -7.812538236389854E-001 -4.590520636167327E-001
             -4.198711506107195E-001 -3.871387741292614E-001
             -2.912685945326537E-001 -2.532555221576769E-001 -1.076727247566867E-001
             -1.025922792354764E+000 -7.812537346931661E-001
             -4.590520549643138E-001 -4.198710925526331E-001 -3.871386455493686E-001
             -2.912685126157339E-001 -2.532553593136832E-001
             -1.076726566570645E-001
                    </eigenvalues>
                    <occupations size="16">
              1.000000000000000E+000  1.000000000000000E+000  1.000000000000000E+000
              1.000000000000000E+000  1.000000000000000E+000
              1.000000000000000E+000  0.000000000000000E+000  0.000000000000000E+000
              1.000000000000000E+000  1.000000000000000E+000
              1.000000000000000E+000  1.000000000000000E+000  1.000000000000000E+000
              1.000000000000000E+000  0.000000000000000E+000
              0.000000000000000E+000
                    </occupations>
                  </ks_energies>
            ```

Finally, comparing the energy of the spin polarized calculation with the spin
degnerate metal calculation, we can see that the spin polarized calculation
gives a lower energy.

!!! pied-piper "Fun facts"
    O2 in its singlet state can be dangerous (see e.g. 
    [`this paper`](https://www.sciencedirect.com/science/article/pii/S1383574211001189)),
    so treating the spin correctly is important!


Iron
----

Now that you've seen how including spin polarization can allow us a correctly
describe the ground state of our system in your calculation, the next step
is to use it to describe a magnetic system.

In a magnetic system there is a net spin polarization in the unit cell. This
means that we'll probably have an odd number of electrons, and the energy of
the system when we include a net spin polarization is lower than the energy
when we don't.

One of the most common magnetic systems is iron, so we'll examine this. The
directory `03_Fe` contains an input file for iron. Note this is a BCC structure
(as set by `ibrav = 3` in the input file), whereas most of the crystals
structures you have examined previously were FCC. The calculation has been set
up in the usual way for a metallic system.

!!! example "Task 3.1 - fixed magnetization"

    1. Run this calculation and check everything worked as expected. What is the
       final energy?

        ??? success "Answer"
            The final energy should be -55.52528610 Ry.

    2. Now make a copy of the calculation directory and in this, modify the
       calculation to turn on spin polarization. Try running the calculation
       with `tot_magnetization = 0.0` first, and compare your total energy to that
       obtained using doubly degenerate bands. 

        !!! note "Note" 
            while in the case of the O2 above, we were able to get our
            calculations to at least converge by using a metallic occupation
            instead of using spin polarization, in the case of iron, it will still
            be a metal when you use spin polarization, so you should not remove
            the input variables associated with this. 

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
the value which gives the lowest overall total energy. However, we can instead
pass an option that tells quantum espresso to automatically find the best
value. This is done by setting the `starting_magnetization` input variable.

!!! example "Task 3.2 - Relaxed magnetization"

    1. Make another copy of the `02_Fe` directory, and this time set `nspin = 2`,
       and `starting_magnetization = 1.0` (do not include the
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
    
    2. See if you can use what we covered in previous labs to calculate and make a
       plot of the electronic band structure of BCC Fe.

        - Plot the spin-up and spin-down bands in different colours.
        - Indicate the Fermi energy on your plot in some sensible way.
        - As the Brillouin zone is different to the ones you have calculated so
          far you'll need to select a few sensible high-symmetry points yourself
          to plot with :slight_smile:.

        ??? success "Answer"
            You can find the relevant input file in the directory 
            `03_Fe/extra_bandstructure`. The band structure should look like The
            following:
            <figure markdown="span">
              ![Diamond primitive cell](assets/Iron_bands.png){ width="500" }
            </figure>

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
