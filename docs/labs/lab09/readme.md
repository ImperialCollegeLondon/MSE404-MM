Charged system and excited states
=================================

**Reminder** Don't forget to copy the `lab09` folder from `/opt/Courses/MSE404/lab09`
to your home directory.

In this lab we're going to see how properties of charged systems can be obtained with DFT,
and how to improve the description of excited states (remember DFT is formally a theory that
allows us to compute ground state and excited states properties by knowing the ground
state density only!).

Ionization potential and electron affinity
------------------------------------------

As a first step, we will look at how to compute ionization energies
and electron affinities in molecules, specifically
carbon monoxide (CO). To do this, we will employ the $\Delta$SCF method, i.e. we will perform
SCF calculations for the neutral molecule (with $N$ electrons) and the charged molecule ($N\pm 1$ electrons)
at the geometry of the neutral molecule and compute the difference in the total energies ($E_{N} \pm E_{N+1}$).
In doing so, we are neglecting the effect of geometry relaxation in the charged system, and therfore the difference in  
nuclear-nuclear repulsion energies, that would be different in two different geometries. However, photoemission is
a fast process and therefore the neglect of relaxations in the charged state is often justified when comparing to 
experiments.
 
Due to the long-range nature of the Coulomb interaction, electrostatic sums on an infinite lattice, or equivalently, 
in a system subject to periodic boundary conditions (PBC), are not always absolutely convergent and in fact 
the sum is divergent if the system is charged, which makes the above-mentioned energy differences ill-defined. 
When using PBC, one can simulate an isolated molecule using a supercell with enough vacuum space. However, the energy calculated for a finite 
supercell differs from that of an infinite supercell, because of the spurious interactions 
of the charged system and its periodic images (as we have seen already this is also true for 
neutral system with no permanent dipole, but in this case the electrostatic sum is absolutely 
convergent and one can get away with a relatively small supercell size).

In DFT code with PBC, charged systems are usually addressed by introducing a compensating 
uniform jellium background, such that the total charge of the supercell is zero. This is of course an approximation 
and it's only justified in the $L\rightarrow\infty$ limit, where $L$ is the linear dimension 
of the supercell. Hence, one way to study charged systems within DFT is to compute the total 
energy for increasing values of the supercell size and extrapolating from these the 
value at $L\rightarrow\infty$. 

More refined methods have been developed to treat charged
systems within PBC and small supercells, such as the Makov-Payne method and Martyna-Tuckerman method, 
which we will also use in this lab. The Makov-Payne method applies a cell size-dependent correction to the total energy
based on an asymptotic analysis of the total energy of charged system with respect to the supercell size. This allows to
get a faster convergence with respect to the supercell size. The Martyna-Tuckerman
method, on the other hand, introduces a cut-off in real space beyond which interactions are set to zero. As a rule of thumb, the supercell 
should be more than twice as large as the size of the molecule for this to be a good approximation.

In Quantum ESPRESSO you can enable the calculation of charged systems by setting
`tot_charge = +1`, if one electron is removed, or `tot_charge = -1`, if one electron
is added, and so on, in the `SYSTEM` section of the input file. When having an odd number
of electrons (e.g. an unpaired electron) we also need to perform a spin-polarised calculation, since
we now have a different number of electrons in the two spin-channels. This can be done by
setting `nspin=2`. Finally, we need to specify the occupation of the Kohn-Sham states. In molecules,
with non-dengerate HOMO-LUMO, we want that each KS state is either fully occupied or empty. This can be
achieved by setting `occupations=fixed`.

The molecular orbitals diagram for the CO molecule can be found at [this website](https://www.chemtube3d.com/orbitalsco/) 

### _Task_

Compute the ionization energy of carbon monoxide. You can find input files and annotated scripts in the
[`01_carbon_monoxide`](01_carbon_monoxide/) folder.

- Inspect the `CO_neutral.in` (neutral) and `CO_charged_p1.in` (positively charged, one electron removed) template files
- Use the `run_cell.sh` bash script to run SCF calculations for increasingly large cell sizes. 
  The script generates a text file `<template_name>_etot_v_box.dat` with two columns: first column is the linear 
  dimension of the cubic cell in angstroms and second column is total energy in Rydberg.
- The script can be used both for the neutral and the charged system. You will have to modify the `template` 
  variable in the script accordingly. Run the script for both neutral and charged system.
- Check it is converged for all cell sizes. 
- After running the script for the neutral and charged molecule, create a text file with two columns, with first 
  column being the linear dimension of the cell in angstroms and second column being the difference 
  between total energies of charged and neutral molecule in Rydberg.
- Use the `deltaE_v_box.gnu` gnuplot file to do a linear fit ($f(x) = a + bx$) of $E_{N-1}-E_{N}$ vs $1/L$. Modify the 
  script such that it reads the file you have generated in the previous point. In the script you will have to substitute
  `FILENAME` with the name of your file.
  Run the script by typing `gnuplot deltaE_v_box.gnu`. 
  This should generate a plot with and  you should get the resulting $a$ and $b$ parameters as well as the standard 
  error from the fitting displayed in the terminal
```bash
Final set of parameters            Asymptotic Standard Error
=======================            ==========================
a               = 14.19            +/- 0.01085      (0.07643%)
b               = -22.8611         +/- 0.1496       (0.6543%)

correlation matrix of the fit parameters:
                a      b      
a               1.000 
b              -0.972  1.000
```
- Extract the value for $L\rightarrow\infty$, i.e. $1/L=0$. The experimental value is $14.01$ eV. How does it compare with
  the DFT result?

### _Task_

Compute the electron affinity of carbon monoxide. 
- Re-do the same steps of previous task replacing `CO_charged_p1.in` with
`CO_charged_m1.in`. 
- For the calculation of the electron affinity we have changed the occupations from `fixed` to `smearing`. 
From the diagram of molecular orbitals of CO, can you explain why?
- The experimental value is $1.326$ eV. How does this compare with the DFT result?

### _Task_
Compute the electron affinity and ionization energy using the Makov-Payne method and Martyna-Tuckerman method to treat 
the charged system at a relatively small cell size. In [`02_carbon_monoxide`](02_carbon_monoxide/) folder copy the input files of the 
charged systems and set the supercell size to 16 Angstroms. The cartesian coordinates for the C and O atoms are
```bash
 C  8.000 8.000 7.436
 O  8.000 8.000 8.564
``` 
To use these methods, we need to add `assume_isolated='mp'` for Makov-Payne and `assume_isolated='mt'` for Martyna-Tuckerman
in the `SYSTEM` section.

Run `pw.x` with these input files and compare the ionization energy and electron affinity with that of 
an infinite supercell from the extrapolation. 

Excited states and band-gap problem - Titanium dioxide (Rutile)
---------------------------------------------------------------
It is well known that in insulators and semiconductors the fundamental band gap is underestimated by DFT with
local (LDA) and semilocal (GGA) exchange-correlation (XC) functionals. This can be related to local and semilocal
XC functionals being unable to remove spurious self-interactions arising from the Hartree term. For instance, if 
we consider a simple one-electron system like the hydrogen atom, one can easily see that the Hartree energy
$E_{H}[n] = \frac{1}{2} \int \frac{n(\mathbf{r})n(\mathbf{r}')}{|\mathbf{r}-\mathbf{r'}|} d^3r d^3r'$ implies
an unphysical self-interaction of the electron with itself. This contribution should be compensated by the
XC term, but an exact cancellation is not possible with local and semilocal functionals. This means that a percentage
of this spurious contribution remains and pushes up the energies of the occupied states. At the same time, unoccupied
states do not contribute to the total density and therefore no self-interaction term arises from them. As a net result,
the gap between occupied and unoccupied states is reduced.

To ameliorate this problem several so-called hybrid XC functionals have been proposed. These are non-local XC functionals
based on the electron density and Kohn-Sham orbitals as well, i.e. $E_{XC}^{hyb} [n(\mathbf{r}),{\phi_{KS}(\mathbf{r})}]$.
The main idea, is to add to the (semi)local XC functionals a percentage of the exact exchange term, which is a functional
of the KS orbitals. In fact, from Hartree-Fock 
theory it is well known that there is a perfect cancellation between the Hartree and the Exchange terms for the self-interaction.

In this part we are going to compute the band structure of titanium dioxide (rutile phase) with a semilocal XC functional (PBE) and 
we will see how using the `B3LYP` hybrid XC functional improves the band gap compared to the experimental value. You are not going to
explicitely run a DFT calculation with a hybrid functional as this is quite computationally demanding, particularly with a plane-wave
basis set (can you think of why?)

The B3LYP energy functional is a very popular hybrid functional and it is obtained by: $E_{XC}^{B3LYP} = (1-a)E_{X}^{LSDA} + a E_{X}^{HF} + b \Delta E_{X}^{B88} + (1-c)E_{C}^{LSDA} + c E_{C}^{LYP}$, where the subscrpits $X$ and $C$ mean exchange and correlation, respectively. The superscripts, on the other hand, identify the different functional forms, for instance $HF$ means Hartree-Fock, $LSDA$ means local spin-density approximation, $B88$ Becke88 functional and $LYP$ the Lee-Young-Parr functional. The three parameters in B3LYP have the following values $a=0.2$, $b=0.72$ and $c=0.81$. Finally $\Delta E$ means a gradient correction to the functional. You can see that in B3LYP one introduces a fraction (0.2) of the Hartree-Fock exchange, which is a non-local functional.

The rutile phase of TiO2 is a direct wide-gap semiconductor, with an experimental bandgap of $3.3$ eV. The unit cell is tetragonal with 
`a=4.58` Ang and `c=2.95` Ang and contains six atoms, two Ti atoms and four O atoms. You can visualise the structure using `xcrysden` as usual.

### _Task_

- Let's first calculate the rutile band structure with a semilocal XC functional (PBE). In [`03_rutile`](03_rutile/) folder you will
find an input file `01_rutile_scf.in` for a SCF calculation and two pseudopotential files `Ti.pbe-sp-van_ak.UPF` and `O.pbe-van_ak.UPF`.
Copy these two files into your `pseudo` folder and run `pw.x` to obtain the ground state density.
- Next, let's compute the band structure on a high-symmetry path by running a non-scf calculation. You will have to use the `01_rutile_nscf.in`
input file.
- Finally, let's plot the bands with `bands.x` as you have done in [Lab 4](https://mse404.gitlab.io/labs/lab04/). You will have to use `01_rutile_bands.in` and plot the result with `gnuplot`. At which high-symmetry point is the direct band gap found?
- Go back to the output of the non-scf calculation and compute the band gap from the lowest unoccupied state and highest occupied states at $\Gamma$.
You should get a band gap of $1.9$ eV, which is $~1/3$ smaller than the experimental value.
- Now look at the output file `02_rutile_scf.in`, which has been obtained with the B3LYP XC functional. At the end of the file you can check the 
difference between the lowest unoccupied state and the highest occupied state. This calculation is not fully converged but you can see that the band gap is now much closer to the experimental value.

-------------------------------------------------------------------------------

Summary
-------

- In this lab we've looked at how to treat charged molecules with DFT and periodic
  boundary conditions. We have also looked at how ionization energies and electron
  affinities can be computed with $\Delta$SCF.
- We've also looked at how more refined methods can help with the convergence
  of the total energy with respect to the cell size for charged molecules, e.g.
  the Makov-Payne method and Martyna-Tuckerman method.
- In the second part we have looked at hybrid exchange-correlation functionals and
  how these can improve the description of excited states and ameliorate the 
  band-bap problem in insulators and semiconductors.

