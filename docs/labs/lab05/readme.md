Forces, Stresses and Structures
===============================
For this week and the next, we will focus on predicting the structural properties of materials. These include finding the stablest structure of molecules and calculating the energy of vibrations in crystals. We will begin by briefly reviewing the very important potential energy surface (PES). 

The potential energy surface
----------------------------
The potential energy surface (PES) is the surface obtained by plotting the total energy of a molecule or crystal against different atomic positions. It is commonly denoted as $U(\{\mathbf{R}\})$, where $\mathbf{R}$ represents the set of atomic positions of all atoms. Many useful quantities associated with a given structure, such as the atomic forces in a molecule, the stress and pressure in a crystal unit cell can be calculated by evaluating the relevant deivatives of the PES.  

Ground state structure of molecules 
-----------------------------------
The ground state structure of a molecule simply refers to the structure with the lowest total energy. This structure is also called the optimal structure or stablest structure of the molecule. Mathematically, this structure is associated with the **global** minimum of the PES. Therefore, the atomic forces in this structure will all vanish.  

Most DFT codes, like Quantum Espresso, moves the atoms around based on the atomic forces until the atomic forces "vanish", i.e. become lower than some cutoff. This process is called **relaxation** or **structural optimisation**. The resulting structure obtained from this procedure is called a relaxed structure. 

In the first part of the lab, we will demonstrate how to calculate the atomic forces and find the optimal structure of a molecule in Quantum Espresso. 

!!! Danger
    When relaxing large molecules, there might be other **local** minima in the PES. These local minima are associated with structures which are higher in the total energy, but the atomic forces also vanish. These structures are called **metastable** structures.

    When finding the stablest structures using DFT, there is always a possibility that your calculation is "trapped" in one of these metastable structures. Whereas this is less likely to happen for small molecules, this can happen with large molecules. You will have to be careful! 

!!! Question
    How can we reduce the risk of mistaking a metastable structure as the ground state structure? 
    ??? success "Answer"
        After obtaining a relaxed structure, add some small and random displacement to the atoms, then redo the relaxation. This will likely bring your molecule to a stabler structure if there is any.  

Forces in molecules
-------------------
The force acting on an atom, $\mathbf{F}$, is defined as the first derivative of the PES, $U$, with respect to the position of that atom, $\mathbf{R}$. Mathematically,
$$
\mathbf{F} = -\nabla_\mathbf{R} U.
$$
It is worth noting that, in order to calculate the atomic forces efficiently, the Hellmann-Feynman Theorem is often invoked in most DFT codes. 

Quantum Espresso calculates the atomic forces in a `pw.x` calculation if you set `tprnfor = .true.` and `tstress = .true.` in the input file. Let's take a closer look at the atomic force calculations in the tasks below. 

In following tasks, we are going to study the atomic forces in a distorted methane molecule, $\mathrm{CH_4}$, where one of the hydrogen atoms (the one sitting on the z-axis above the $\mathrm{C}$ atom) is pushed closer to the carbon atom than others. 

!!! example "Task 1 - Examining atomic forces in the output file"
    - Read and run the input file `01_methane_force.in`. Note how the two settings `tprnfor` and `tstress` have been set to true. 

    - When the output file is out, search for the line saying `Forces acting on atoms`. This should be written just before the timing information. 

    - You should find a section that looks like the following: 
    ```python 
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

    - Why are these atomic forces expected for this distorted methane molecule? 
    ??? success "Answer"
        Since the $\mathrm{H}$ atom (atom 2) above the $\mathrm{C}$ atom (atom 1) is approaching the $\mathrm{C}$ atom and the other $\mathrm{H}$ atoms below the $\mathrm{C}$ atom (atoms 3 to 5), the upper $\mathrm{H}$ atom should be repelled by the rest of the atoms. Therefore, atom 2 should experience an upward force and the rest should experience a downward force. 

        This is indeed the case. The force on atom 2 is vertically upwards, and the forces on the rest of the atoms are vertically downwards (the atomic forces along the $x,y$-directions are negligible). You can also check that all the $z$-components of the downward forces acting on atoms 1 and 3 to 5 sum up to the same magnitude as the $z$-component of the upward force acting on atom 2.  

    - What are the units of atomic forces in Quantum Espresso? 
    ??? success "Answer"
        Ry/Bohr. 1 Ry/Bohr is equivalent to ... eV/Å (check this yourself!). In this task, the $z$-component of the force acting on the $\mathrm{H}$ atom is 0.0522 Ry/Bohr. As a rough approximation, this means drifting the $\mathrm{H}$ atom in the positive $z$-direction by 1 Bohr should lead to a reduction in total energy of about 0.0522 Rydberg. This is equal to 0.710 eV, which is quite a lot.

!!! Note 
    The `Total force` listed here is the square root of the sum of _all_ of the force components squared rather than the sum of the magnitudes of the individual forces on the atoms. I'm not sure why, but it's likely because the number is intended more as a guide to check overall accuracy. If the `Total SCF correction` is comparable to the `Total force` it usually means you need to try to better converge the SCF cycle (via `conv_thr`).


Just as we have to check the convergence of the total energy against the PW cutoff, we have to check the convergence of atomic forces against the PW cutoff too. Let's do this in the next task.  

!!! example "Task 2 - Examining convergence"
    In this task, multiple input files of the same distorted methane molecule are prepared. Each input file is labelled by their PW cutoff. 

    - Run the `whatever auto.sh` script to run single-point calculations for all input files. 

    - Run the `whatever python or bash` to obtain the convergence curves. 

    - Which quantities are being plotted by the script for determining the convergence? 

    ??? success "Answer" 
        The script plots the fractional difference of the total energy with respect to the best-converged total energy, $(E_{\mathrm{best}}-E_{i})/E_{\mathrm{best}}$, and analogously the fractional difference in the total force $(F^{\mathrm{tot}}_{\mathrm{best}}-F^{\mathrm{tot}}_{i})/F^{\mathrm{tot}}_{\mathrm{best}}$.   

Optimisation of molecular geometry 
-------------------
In addition to calculating the atomic forces associated with a particular structure, another important functionality of Quantum Espresso is optimising the structure of molecules, which we will try out in this following task.

!!! example "Task 3 - Optimising the molecular structure"
    - You'll need to add an `IONS` section to the input. This is the only mandatory addition. There are several variables that can be specified within this section (it can be empty if you're happy with all defaults), which control the algorithm used to find the optimal ionic positions. Consult the PW documentation for details.


!!! note "Interlude"
    A short break you should take, tired if you are. 


Ground state structure of crystals 
-----------------------------------
We have seen how the stablest structures of molecules can be found by minimising the atomic forces in DFT. In principle, it is straightforward to apply the idea for the optimisation of the structures of crystals. However, an additional consideration arises from optimising the lattice constant of the crystal alongside the atomic positions within the unit cell. This leads to the need for a few more quantities 

Stress in crystals
-------------------
Stress, $\sigma_{ij}$, can be intuitively understood as the ratio between the $i$-th component of the induced force (e.g. this could be the $x$-component), $F_{i}$ , due to a deformation along the $j$-th direction (e.g. in the $x$- or $z$-direction), $\epsilon_{j}$. Therefore, mathematically, the stress is a tensor and is defined through the relation  

$$
F_{i} =  \sigma_{ij}\epsilon_{j}
$$

We will see later below how the deformation of a crystal along one direction induces forces along other directions. 

Pressure and bulk modulus of crystals 
-------------------
Whereas stress describes the force induced on a crystal due to a distortion along a certain direction, the pressure and the bulk modulus of a crystal measures the changes in total energy due to contractions or expansions. Pressure, P, is defined as the negative of the first derivative of the PES against the volume, V. It measures the tendency of the crystal to expand or contract. Mathematically,  

$$
P = -\frac{dU}{dV}.
$$

If the pressure is negative, the crystal structure is unstable and will contract. If the pressure is positive, the crystal structure is also unstable but it will expand.   

The bulk modulus per volume is defined as the negative of the second derivative of the PES against the volume. It measures the stability of the crystal against expansion or contraction. Mathematically,  

$$
B = -V\frac{d^2U}{dV^2}.
$$

If the bulk modulus is negative, the curvature of the PES against volume is positive. This means the structure is stable and it costs energy to compress or expand the crystal. Having introduced these useful quantities, we can take a look at these quantities in the output file of Quantum Espresso. 

In the following tasks, we will look at the poly(para-phenylene) (PPP) polymer. It is a chain of benzene rings which extends essentially infinitely in one direction. It is therefore a one-dimensional crystal. The unit cell of this crystal contains two benzene rings. We will first run a singe-point calculation of this crystal.  


!!! example "Task 4 - calculating stress and pressure of crystals"
    
    - Read the input file `04_PPP.in`. 
        ```python 
        
        K-GRID LINE   #(1)! 
        
        ```
    
        1.   Note that here we have added the k-grid flag, which is required for a calculation of crystals. We have chosen the optimal k-grid for you.  
        
    - Look for the line saying "stress"

    - What are the units for the stress tensor?
    ??? success "Answer"
        Ry/Bohr**3. The dimension of this unit is energy per volume, which is the same dimension as pressure. 

    - What are the units for pressure? 
    ??? success "Answer"
        kPa. This is one of the common units for pressure. 

If you take a closer look at the PPP molecule again, you will notice a torsion angle, $\theta$, between two neighoburing rings in the unit cell. Our next task is to optimise this torsion angle.  

Optimisation of the unit cell of crystals 
-------------------
The first method of finding the optimal torsion angle is to find the minimum in the PES, plotted against the torsion angle $\theta$. This is what we will do in the following task. 

!!! example "Task 5 - Optimising structure through PES"
    - Multiple input files have been prepared for different torsion angles. A bash script called "auto_torsion.sh" have been prepared for you to run. This script automatically runs all the `pw.x` for all the input files. 

    - Once all the input files have been run, run the Python script `theta.py`, which plots the total energy of the PPP molecule against the torsion angle. Read the Python script and make sure you understand what it is doing. 

    ??? success "Result"
        A plot will be uploaded very soon 

The second method is to just run a "relax" calculation with Quantum Espresso. 
!!! example "Task 6 - Optimising structure through DFT relaxation"
    - Go to directory 06_PPP

    - Open the input file named 06_PPP. The torsion angle of the initial structure in this input file is 20 degrees. Edit the relevant line so that this calculation is a relaxation calculation. 
    ??? success "Result"
        ` calculation = relax `
        `&IONS
        /`


    - Open the output file in `xcrysden`. As you will see, it's possible to just open the final optimized structure, or load the optimization as an animation. Since we already started quite close to the final structure, you won't see much of a change if you watch an animation. 

    - Use the `dihedral angle` option in xcrysden to find the torsion angle of the relaxed   structure. How does this compare to what you previously predicted?


Optimizing Unit Cells
---------------------


In doing this you should keep in mind that it can take a higher energy-cut off to converge the stress as we saw at the start of this lab. Also, if you recall the number of planewaves depends on the cell volume, so if during the course of the calculation we change the volume we may also change the number of plane waves, or if we fix the number of plane waves we are changing the effective energy cut-off (the latter for Quantum Espresso, but different codes will handle this differently). So you need to be quite careful about this when optimizing lattice vectors.

In Quantum Espresso we can do this variable-cell relaxation by setting `calculation = 'vc-relax'` in the `CONTROL` section. We must additionally specify both an `IONS` section as previously, along with a `CELL` section.

!!! example "Task 7 - Optimising the unit cell"
    - Read the input file 07_PPP. You'll notice in addition to the inputs mentioned, there's also a fairly high energy cut-off, and we've lowered the SCF convergence threshold from the default. Try running this now and let's look at the output.

    !!! warning 
        The output is a little different in this case, since at the end of the optimization, an `scf` calculation is automatically performed starting from the optimized structure. This is because Quantum Espresso fixes the basis set as that for the original input structure during the calculation, so if the structure has changed a lot, you may calculate a different stress when starting from the relaxed structure. You will be able to see from the final stress or pressure whether you should rerun your calculation.

!!! note 
    Additionally, the cell output will likely be in Quantum Espresso's representation where the cell vectors are given explicitly in Bohr along with a scaling factor `alat` which is fixed. (In our case here `alat` will be `A` converted from Angstrom to Bohr). If you want to rerun your calculation you could either input the cell using these directly, or calculate appropriate values for the input. You may need to do the latter if you want to find the new lattice length anyway.

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

