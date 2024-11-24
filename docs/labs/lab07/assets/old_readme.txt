Forces, Stresses and Structures
===============================

The contributor to this tutorial is currently getting used to mkdocs. So the following are the tentative structure of the updated tutorial and serve as practices for the contributor. 

!!! danger 
    I am currently figuring out the best way to modify the workflows of each task. 

!!! warning 
    mkdocs is amazing. Thanks Chengcheng.  
------------------------- Tentative actual text ---------------------

For this week and the next, we will focus on predicting the structural properties of materials. These include the stablest geometry of molecules and the energy of phonons in crystals. We will begin by briefly reviewing how the potential energy surface (PES) is very useful for these tasks. 

The potential energy surface
----------------------------


Forces and stress in molecules
-------------------


!!! example "Task 1 - Examining convergence"
    The details of the task are to be modified  

    This is the second line of the task 

    - Will we make some sort of plot of total force against KE? 

        ??? success "Answer"
            Yes, I think we really should.  

    - Will a new python script replace the bash script for energy convergence? 
        
        ??? success "Answer"
            Yes. I will make that asap.

!!! warning 
    I think it will be helpful to remind readers to be careful of the use of different units, fractional coordinates, etc. in Quantum Espresso


!!! note 
    For job submission, I think it is the easiest to use auto.sh. 

    For data analysis or output file reading/input file modifications, it is the easiest to use Python. 


Optimisation of molecular geometry 
-------------------

Constrained relaxation 
-------------------

!!! Quiz 
    
    - What are the advantages of constraining the motion of some atoms?
        
        ??? success "Answer"
            
            The advantages includes reducing the number of relaxation steps required for reaching the equilibrium structure. 


Pressure and bulk modulus of crystals 
-------------------

Optimisation of the unit cell of crystals 
-------------------


------------------------- old content below -------------------------

**Reminder** Don't forget to copy the `lab05` folder from `/opt/Courses/MSE404/lab05`

To find the minimum energy position of an atom, we could manually move it,
calculating the total energy each time, effectively finding the position where
the total force on it is zero. Some examples of how you can do this are given
[here](../extras/labs/using_total_energies/readme.md), you can take a look
through in your own time if you're interested. Instead, in this lab we'll be
looking at how forces and stresses can be calculated within DFT at essentially
no extra cost via the Hellman-Feynman theorem, and how these can be used.

In Quantum Espresso you can enable the calculation of forces and stresses by
setting `tprnfor = .true.` and `tstress = .true.` respectively in the
`CONTROL` section of the input file.

Forces in Methane
-----------------

As a first example, let's look at methane and calculate how the forces
converge with planewave energy cut-off. We have set up a template input file
as before in
[`01_forces/01_methane/CH4_base.in`](01_forces/01_methane/CH4_base.in). Take a
look at this now. The only new settings here are the two additional variables
mentioned above.

We also have a simple script to run this template file with energy cut-off
values from 20 to 60 Ry as [`auto_run.sh`](01_forces/01_methane/auto_run.sh).
Run this now and take a look at one of the output files. You'll see before the
final timing information a section that looks like the following:

```python
     Forces acting on atoms (cartesian axes, Ry/au):

     atom    1 type  1   force =     0.00000620    0.00000000    0.00002841 
     atom    2 type  2   force =     0.00000137    0.00000000    0.01218625
     atom    3 type  2   force =     0.01151561    0.00000000   -0.00407450
     atom    4 type  2   force =    -0.00576159   -0.00997884   -0.00407008
     atom    5 type  2   force =    -0.00576159    0.00997884   -0.00407008

     Total force =     0.024421     Total SCF correction =     0.000247     #(1)!


     Computing stress (Cartesian axis) and pressure

          total   stress  (Ry/bohr**3)                   (kbar)     P=   -0.06
  -0.00000044   0.00000000   0.00000000         -0.06      0.00      0.00
   0.00000000  -0.00000044   0.00000000          0.00     -0.06      0.00
   0.00000000   0.00000000  -0.00000044          0.00      0.00     -0.06
```

1. Here the square root of the sum of the squares of each component of the forces on each atom is being shown. See below for more explanation 

So we have a list of the components of the force acting on each atom along
Cartesian axes, then a total force and total scf correction.
**Note** the `Total force` listed here is the square root of the sum of _all_
of the force components squared rather than the sum of the magnitudes of the
individual forces on the atoms. I'm not sure why, but it's likely because the
number is intended more as a guide to check overall accuracy. If the
`Total SCF correction` is comparable to the `Total force` it usually means
you need to try to better converge the SCF cycle (via `conv_thr`).

Following this we can see the calculated stress and corresponding pressure
on the unit cell. While this number doesn't mean much in principle for a
calculation on a molecule in a box, if these numbers are not all close to
zero it indicates your box size is likely not big enough.

### Convergence

Now let's look at the convergence of the forces. As we've been doing in
previous labs, we can extract the total force as a function of the energy
cut-off using [`awk`](../extras/misc/linuxcommands/readme.md#awk):

```bash
awk '/kinetic-energy/{ecut=$4} /Total force/{print ecut, $4}' *out
```
This works by saving (space delimited) field number 4 as a variable `ecut` on
a line containing the string `kinetic-energy`, and then on a line containing
the text `Total force` it outputs the value of this variable along with the
text in the fourth field.

We could modify this to also output the total energy as
```bash
awk '/kinetic-energy/{ecut=$4}
     /!.*total/{etot=$5}
     /Total force/{print ecut, etot, $4}' *out
```
(You can add line breaks for clarity within an awk command if it gets long).

And add in the total pressure also with
```bash
awk '/kinetic-energy/{ecut=$4}
     /!.*total/{etot=$5}
     /Total force/{totfor=$4}
     /total.*stress/{print ecut, etot, totfor, $6}' *out
```

Try running this and saving the convergence to a file.

### _Task_

- Plot the fractional difference with respect to the most well converged
  result, i.e. $|\frac{x_{conv} - x_i}{x_{conv}}|$, for each of total energy, 
  total force and pressure as a function of energy cut-off.
- What can you see about the rate of convergence of each of these parameters?


Optimizing Ionic Positions - PPP
--------------------------------

Now let's look at a polymer, in this case poly(para-phenylene) (PPP), which consists
of a chain of benzene rings, with a torsion angle of $\theta$, which should be around
$30^{\circ}$. If we wanted to predict the value of $\theta$, we could take structures
with a few different angles and find the minimum of the energy by fitting a parabola. 
In the directory `02_structure/01_ppp` there are some input files for different torsion
angles, where e.g. [`02_structure/01_ppp/PPP_30.in`](02_structure/01_ppp/PPP_30.in) corresponds to
$\theta=30^{\circ}$. Since the polymer is periodic in only 1 dimension (the z axis),
we have added empty space in the other two directions. We have also added k-point sampling
along the z-axis only.

### _Task_
- Run the input files for PPP for angles of 20, 25, 30, 35 and 40 degrees. Try writing a script
  to automate the process. Your script should generate an output file called `e_v_theta.txt`
  which has two columns - the first should be the torsion angle and the second should be the
  calculated energy.
- Use the file `plot_ppp.gp` to plot energy vs. $\theta$. This will generate a file called
  `ppp.eps`, where a quadratic function has been fit to the data. Using this fit, work out
  what the optimum torsion angle is.

While the above approach gives us a good idea of the optimum torsion angle, it's a bit tedious
to generate the input structures, and it doesn't tell us whether or not the benzenes remain
rigid, or if there are some internal distortions. To tell us this, we can use the forces.
In fact, given that we can readily calculate forces, wouldn't it be nice if the code could
automatically use these forces and find the atomic positions where these forces are zero
(or at least within some tolerance of zero)?

In Quantum Espresso this type of calculation, where the atomic positions are
relaxed to their minimum energy positions can be performed by selecting
`calculation = 'relax'` in the `CONTROL` section. There are a number of
additional variables that you can specify to control this process:

- `tprnfor` is automatically set to `.true.` so that forces are computed.
- You'll need to add an `IONS` section to the input. This is the only
  mandatory addition. There are several variables that can be specified within
  this section (it can be empty if you're happy with all defaults), which
  control the algorithm used to find the optimal ionic positions. Consult the
  PW documentation for details.

### _Task_

The directory `02_structure/01_ppp` also contains an input file
[02_structure/01_ppp/PPP_opt.in](02_structure/01_ppp/PPP_opt.in) for relaxing the
structure of PPP, where we start from the $\theta=30^{\circ}$ structure since this
is close to the minimum. Run this and take a look at the output file once the 
calculation finishes - it might take a couple of minutes. The interesting bit
is right near the end, just before the timing output. It should look something
like the following:
```
     Forces acting on atoms (cartesian axes, Ry/au):

     atom    1 type  1   force =     0.00001266    0.00000665   -0.00024259
     atom    2 type  1   force =     0.00066443    0.00046526    0.00036559
     atom    3 type  1   force =    -0.00064066   -0.00045488    0.00039018
     atom    4 type  1   force =     0.00066443    0.00046526   -0.00036559
     ....
     atom   17 type  2   force =     0.00017760    0.00095599   -0.00028356
     atom   18 type  2   force =    -0.00020428   -0.00097249   -0.00028734
     atom   19 type  2   force =     0.00017760    0.00095599    0.00028356
     atom   20 type  2   force =    -0.00020428   -0.00097249    0.00028734

     Total force =     0.003777     Total SCF correction =     0.000132

     bfgs converged in   6 scf cycles and   5 bfgs steps
     (criteria: energy <  1.0E-04 Ry, force <  1.0E-03 Ry/Bohr)

     End of BFGS Geometry Optimization

     Final energy   =    -144.7886351866 Ry
Begin final coordinates

ATOMIC_POSITIONS (bohr)
C        6.140945672   6.140862654   1.381651769
C        6.169063980   3.869814143   2.732426642
C        6.112760385   8.411811375   2.732410276
C        6.169063980   3.869814143   5.337957358
....
H        8.194183393   2.617484701   9.770830680
H        4.087214702   9.663976947   9.770881223
H        8.194183393   2.617484701  14.440322320
H        4.087214702   9.663976947  14.440271777
End final coordinates
```

- We have the force output. They should all be pretty close to zero now.
- Then it should say that our relaxation converged (`bfgs` is the name of the
  default algorithm used) and how many steps it took.
- We have the final total energy output.
- Finally we have the optimized atomic positions. 

In order to be sure we have accurately relaxed the structure, we should converge
the forces with respect to cut-off energy, k-points and x and y cell dimensions.
There are also some variable which affect when the calculation stops. These include

- `etot_conv_thr` and `forc_conv_thr`, which are used to determine when the optimization
  finishes (listed following `criteria:` in the output above) and
- `conv_thr`, which as we've already seen is the scf convergence criteria.

### _Task_

- Open the output file in `xcrysden`. As you will see, it's possible to just open the
  final optimized structure, or load the optimization as an animation. Since we already
  started quite close to the final structure, you won't see much of a change if you
  watch an animation.
- Use the `dihedral angle` option in xcrysden to find the torsion angle of the relaxed
  structure. How does this compare to what you previously predicted?


Fixing Some Atoms - Methane
---------------------------

Finally, it's useful to know that we can also optimize the positions of selected atoms
within a system, keeping some atoms fixed.  The input file
[02_structure/02_methane/CH4.in](02_structure/02_methane/CH4.in) shows how this works
for the example of methane.

### _Optional Task_

- Take a look at the input file.
- We've added some 1s and 0s following the atomic positions. These define
  multiplying factors for the various calculated force components on an atom.
    - By adding three 0s following the carbon position we ensure it will be
      fixed at 0,0,0.
    - We only allow the first H atom to move along the z-axis.
    - And we only allow the second H atom to move within the x-z plane.
- Run `pw.x` with this input file and check the result.
- Has the C-H bond been lengthened or shortened?
- How much lower is the total energy compared to the starting configuration?
- Perform this calculation for energy cut-offs of 10 and 50 Ry, and see how
  this affects predicted C-H bond length.
- What do you think would happen if you did this calculation for carbon
  diamond? What about graphite?


Optimizing Unit Cells
---------------------

In a similar way to using the calculated forces to optimize the ionic
positions we can also use the calculated stresses to optimize the unit cell.
In doing this you should keep in mind that it can take a higher energy-cut off
to converge the stress as we saw at the start of this lab. Also, if you recall
the number of planewaves depends on the cell volume, so if during the course
of the calculation we change the volume we may also change the number of plane
waves, or if we fix the number of plane waves we are changing the effective
energy cut-off (the latter for Quantum Espresso, but different codes will
handle this differently). So you need to be quite careful about this when
optimizing lattice vectors.

In Quantum Espresso we can do this variable-cell relaxation by setting
`calculation = 'vc-relax'` in the `CONTROL` section. We must additionally
specify both an `IONS` section as previously, along with a `CELL` section.

An example [input file](02_structure/03_silicon/Si.in) for silicon is given in the
directory [`02_structure/03_silicon`](02_forces/03_silicon). Take a look at this
now. You'll notice in addition to the inputs mentioned, there's also a fairly
high energy cut-off, and we've lowered the SCF convergence threshold from
the default. Try running this now and let's look at the output.

The output is a little different in this case, since at the end of the
optimization, an `scf` calculation is automatically performed starting from the
optimized structure. This is because Quantum Espresso fixes the basis set as
that for the original input structure during the calculation, so if the
structure has changed a lot, you may calculate a different stress when
starting from the relaxed structure. You will be able to see from the final
stress or pressure whether you should rerun your calculation.

Additionally, the cell output will likely be in Quantum Espresso's
representation where the cell vectors are given explicitly in Bohr along with
a scaling factor `alat` which is fixed. (In our case here `alat` will be `A`
converted from Angstrom to Bohr). If you want to rerun your calculation you
could either input the cell using these directly, or calculate appropriate
values for the input. You may need to do the latter if you want to find the
new lattice length anyway.

### _Task_

- So for the silicon case, you'll see the components of the lattice vectors
  have changed from 0.5 to 0.497277 or something close to that. What is our
  predicted silicon lattice constant?

------------------------------------------------------------------------------

Summary
-------

In this lab we have seen

- How to output forces and stresses in our calculation, and how to check
  these quantities have converged.
- How to use these forces in a calculation to optimize the atomic positions,
  where they are moved until the forces on the atoms are less then some
  threshold value.
- How the stresses can be used to optimize the structure, where the lattice
  constant is changed until the stress on the system is less than some
  threshold value.

------------------------------------------------------------------------------


### Extra: Optimizing Graphene

The folder [`02_structure/04_graphene`](02_structure/04_graphene) contains an
[input file](02_structure/04_graphene/C_graphene.in) for graphene. The thing
to note in this case is the additional use of `cell_dofree = '2Dxy'` in the
`CELL` section. We have used this to say that we only want to optimize the
cell in the xy plane. The spacing between periodic images in the z-direction
should be large enough such that the interaction is small, but we otherwise
don't care about the stresses in this direction.

### _Task_

- Starting from the provided graphene input file, find the length of the C-C
  bond in graphene to 3 significant figures.
    - Be sure to test your result is converged with respect to plane wave
      energy cut-off and k-point sampling, along with the internal tolerances
      being sufficiently small.

