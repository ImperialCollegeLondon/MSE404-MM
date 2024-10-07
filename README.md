Modelling Materials with Density Functional Theory
==================================================

## Introduction
This is the computer lab material for the course "Modelling Materials with Density Functional
Theory" taught in the Department of Materials at Imperial College London.

This course will introduce students to the modelling of materials
with density-functional theory (DFT). In the labs we will use the free,
open-source DFT code [Quantum Espresso](http://www.quantum-espresso.org/). While the format of the input files may be differeint in other DFT codes, the general
principles will be the same. The labs will also briefly introduce students to
the use of the Linux operating system and how it enables the effective use of computational
resources, as well as showing students some simple scripting, analysis, and
visualization tools.

<!-- The labs are set up assuming students are using our student server remotely -->
<!-- which has all the necessary software installed already. If it is not possible -->
<!-- for you to do this, I suggest installing xubuntu to a virtual machine on your -->
<!-- laptop. If you need to do this, there are some guidelines at -->
<!-- [vmsetup](labs/extras/misc/vmsetup/readme.md).  -->

<!-- Note that the remote desktop -->
<!-- software we use, [x2go](https://wiki.x2go.org) is freely available for Windows, -->
<!-- Mac, and Linux so if you'd prefer to use your own laptop, please go ahead. -->

## Structure of the Labs

The labs are organised as follows:

- [Lab 1: Getting started in Linux](docs/labs/lab01/readme.md)
- [Lab 2: Quantum Espresso Input and Output](docs/labs/lab02/readme.md)
- [Lab 3: Converging your Calculations](docs/labs/lab03/readme.md)
- [Lab 4: The Electronic Band Structure and Density of States](docs/labs/lab04/readme.md)
- [Lab 5: Metals and Density of States](docs/labs/lab05/readme.md)
- [Lab 6: Spin Polarization and Magnetic Systems](docs/labs/lab06/readme.md)
- [Lab 7: Structural Optimisation](docs/labs/lab07/readme.md)
- [Lab 8: Vibrational Normal Modes and Phonons](docs/labs/lab08/readme.md)
- [Lab 9: Finite Temperature Properties](docs/labs/lab09/readme.md)

Additional Material:

- [Calculating Useful Properties from Total
  Energies](docs/labs/extras/labs/using_total_energies/readme.md)
- [Examining the Charge Density and Wavefunction](docs/labs/extras/labs/visualising_output/readme.md)
- [The Hydrogen Atom and Electron Spin](docs/labs/extras/labs/hydrogen_atom/readme.md)
- [The Tetrahedron Method for Densities of States](docs/labs/extras/labs/tetrahedron_method/readme.md)
- [Optical Properties and Time-Dependent Density Functional Theory](docs/labs/extras/labs/tddft/readme.md)
- [Charged system and excited states](docs/labs/extras/labs/charged_system/readme.md) 

Extras:

- [Setting up a virtual machine](docs/labs/extras/misc/vmsetup/readme.md)
- [More Useful Linux Commands](docs/labs/extras/misc/linuxcommands/readme.md)
- [Plotting with Gnuplot](docs/labs/extras/misc/gnuplot/readme.md)
- [Bash and Shell Scripting](docs/labs/extras/misc/shellscripting/readme.md)
- [Running in Parallel](docs/labs/extras/labs/running_in_parallel/readme.md)

## How to Use this Material

The latest version of the course is available online at
<https://imperialcollegelondon.github.io/MSE404-MM/>. The text source files named
`readme.md` in the various directories use markdown. You can download the
whole repository and view them in a terminal, or read them at the Github site
<https://imperialcollegelondon.github.io/MSE404-MM/>.


Each lab is designed to be self-contained. The webpage will show the general
direction and explain the concepts, while the input files and scripts are
stored in directories located along side the `readme.md` file, labeled by their
apparences in the text. 

There will be code blocks like the one below in the text
```bash
echo "Hello World"
```
which are meant to be run in the terminal. And you can click the 
:material-content-copy: icon on the right to copy the code to your clipboard.

There will also be inline annotations (⊕) like the one 
below:

```fortran
program
    print *, "Hello World" ⊕
end program hello
```
If you click the ⊕ button, a pop up window will show up (this only works in the 
web page).

There will also be admonitions like the one below:

> [!WARNING] 
> Don'f forget to add comments to your code to explain what you are doing!

where important/optional information is given. Also, the Tasks will also be
marked as admonitions like the one below:

> [!Note] 
> Read this `readme.md` file.

## Admin resources

- [Contributing](admin/contributing.md)
- [Server config](admin/server_config.md)

## Acknowledgements
- Original materials provided by Éamonn Murray (https://gitlab.com/eamonnmurray/MaterialsModelling)
- Contributors: Simão João, Christopher Cheung, Jordan Edwards, Chengcheng Xiao, Indrajit Maity, Valerio Vitale, Laura Ratcliff and Johannes Lischner.
- Webpage refreshed by Chengcheng Xiao, 2024.

