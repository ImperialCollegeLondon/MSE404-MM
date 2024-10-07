Structural Optimisation
===============================
For this week and the next, we will focus on predicting the structural properties of materials. These include finding the stablest structure of molecules and crystals, and calculating the vibrational patterns in molecules and crystals.

In this week, we will learn how to calculate the key quantities for describing the stability of structures of moleculs and crystals. You will have learnt them from the lectures, and we will give a brief review for them in this lab before the calculations to refresh your memory. We will begin by briefly reviewing the very important potential energy surface (PES). 

## The potential energy surface
The potential energy surface (PES) is the surface obtained by plotting the total energy of a molecule or crystal against different atomic positions. It is commonly denoted as $U(\{\mathbf{R}\})$, where $\mathbf{R}$ represents the set of atomic positions of all atoms. Many useful quantities associated with a given structure, such as the atomic forces in a molecule, the stress and pressure in a crystal unit cell can be calculated by evaluating the relevant deivatives of the PES.  

## Ground state structure of molecules 
The ground state structure of a molecule simply refers to the structure with the lowest total energy. This structure is also called the optimal, equilibrium, or stablest structure of the molecule. Mathematically, this structure is associated with the **global** minimum of the PES. Therefore, the atomic forces in this structure will all vanish.  

Most DFT codes, like Quantum Espresso, moves the atoms around based on the atomic forces until the atomic forces "vanish", i.e. become lower than some cutoff. This process is called **relaxation** or **structural optimisation**. The resulting structure obtained from this procedure is called a relaxed structure. 

In the first part of the lab, we will demonstrate how to calculate the atomic forces and find the optimal structure of a molecule in Quantum Espresso. 

!!! danger "Danger: metastable structures"
    When relaxing large molecules, there might be other **local** minima in the PES. These local minima are associated with structures which are higher in the total energy, but the atomic forces also vanish. These structures are called **metastable** structures.

    When finding the stablest structures using DFT, there is always a possibility that your calculation is "trapped" in one of these metastable structures. Whereas this is less likely to happen for small molecules, this can happen with large molecules. You will have to be careful! 

!!! Question
    How can we reduce the risk of mistaking a metastable structure as the ground state structure? 
    ??? success "Answer"
        After obtaining a relaxed structure, add some small and random displacement to the atoms, then redo the relaxation. This will likely bring your molecule to a stabler structure if there is any.  

### Forces in molecules
The force acting on an atom, $\mathbf{F}$, is defined as the first derivative of the PES, $U$, with respect to the position of that atom, $\mathbf{R}$. Mathematically,
$$
\mathbf{F} = -\nabla_\mathbf{R} U.
$$
It is worth noting that, in order to calculate the atomic forces efficiently, the Hellmann-Feynman Theorem is often invoked in most DFT codes. 

In following tasks, we will go through how the atomic forces are can be calculated using Quantum Espresso. We have prepared a distorted methane molecule, $\mathrm{CH_4}$, where one of the hydrogen atoms (the one sitting on the z-axis above the $\mathrm{C}$ atom) is pushed closer to the carbon atom than others. 

!!! example "Task 1 - Examining atomic forces in the output file"
    - Read input file `CH4.in`. The most important section of the file is shown below. 
    ``` python 
     &CONTROL
        pseudo_dir = '.'  
        disk_io = 'none'  
        tprnfor = .true.   #(1)!
     /
    ...
    
    K_POINTS gamma         #(2)!

    ```
        1. `tprnfor` asks Quantum Espresso to print the forces.
        2.  We have to carry out $\Gamma$-point calculations for a molecule. 
    - Run the command `pw.x < CH4.in > CH4.out`.
    - When the output file is out, search for the line saying `Forces acting on atoms`. This should be written just before the timing information. 

    - You should find a section that looks like the following: 
    ```python 
     Forces acting on atoms (cartesian axes, Ry/au):

     atom    1 type  1   force =    -0.00000802    0.00000000   -0.05271243
     atom    2 type  2   force =     0.00000254    0.00000000    0.05256575
     atom    3 type  2   force =    -0.00352803    0.00000000    0.00004914
     atom    4 type  2   force =     0.00176676    0.00306145    0.00004877
     atom    5 type  2   force =     0.00176676   -0.00306145    0.00004877

     Total force =     0.074694      Total SCF correction =     0.000393
    ```

    - Why are these atomic forces expected for this distorted methane molecule? 
    ??? success "Answer"
        Since the $\mathrm{H}$ atom (atom 2) above the $\mathrm{C}$ atom (atom 1) is approaching the $\mathrm{C}$ atom, the covalent $\mathrm{C}-\mathrm{H}$ bond between these two atoms is compressed. The two atoms will repel each other to uncompress the bond. Atom 2 should experience an upward force and atom 1 should experience a downward force. This is indeed the case in the output file.  


    - What are the units of atomic forces in Quantum Espresso? 
    ??? success "Answer"
        Ry/Bohr. (1 Ry/Bohr is equivalent to 25.7 eV/Å, check this yourself!). In this task, the $z$-component of the force acting on the $\mathrm{H}$ atom is 0.0527 Ry/Bohr. As a rough approximation, this means drifting the $\mathrm{H}$ atom in the positive $z$-direction by 1 Bohr should lead to a reduction in total energy of about 0.0527 Rydberg. This is equal to 0.710 eV, which is quite a lot. This is because the covalent bond is quite compressed. 

!!! Note 
    The `Total force` listed here is the square root of the sum of _all_ of the force components squared rather than the sum of the magnitudes of the individual forces on the atoms. It is likely because the number is intended more as a guide to check overall accuracy. If the `Total SCF correction` is comparable to the `Total force` it usually means you need to try to better converge the SCF cycle (via `conv_thr`).

### Convergence test 
Just as we have to check the convergence of the total energy against the PW cutoff, we have to check the convergence of atomic forces against the PW cutoff too. Let's do this in the next task.  

!!! example "Task 2 - Examining convergence"
    In this task, multiple input files of the same distorted methane molecule are prepared. Each input file is labelled by their PW cutoff. 

    - Run `python3 ecutwfc_file_build.py` script to produce a list of input files. This script is very similar to the one used in Lab 3. This will create input files with PW cutoffs from 10 to 100 Ry.  

    - Run `bash ecutwfc_force.sh` to run all the `pw.x` calculations as in Lab 3. 

    - Which quantities are being plotted by the script for determining the convergence? 

    ??? success "Answer" 
        The script plots the fractional difference of the total energy with respect to the best-converged total energy, $(E_{\mathrm{best}}-E_{i})/E_{\mathrm{best}}$, and analogously the fractional difference in the total force $(F^{\mathrm{tot}}_{\mathrm{best}}-F^{\mathrm{tot}}_{i})/F^{\mathrm{tot}}_{\mathrm{best}}$.   
    
    - How does the convergence curves look like? Which quantity converges faster with PW cutoff? 
    ??? success "Answer"

### Optimisation of molecular geometry 
Now that we have seen how to calculate the atomic forces associated with a particular structure, we can use Quantum Espresso to optimise the structure of molecules.

!!! example "Task 3 - Optimising the molecular structure"
    - Read the input file `CH4_opt.in`, the important parts are shown below. 
    ```python
     &CONTROL
        calculation = 'relax'   #(1)! 
        pseudo_dir = '.'
        disk_io = 'none'
        tprnfor = .true.
        forc_conv_thr = 1.0D-4    #(2)!
     /
    
     &SYSTEM
        ibrav =  1
        A = 20.0
        nat = 5
        ntyp = 2
        ecutwfc = 80
     /
    
     &ELECTRONS
     /
    
     &IONS        #(3)!
     / 
    ...
    ```
        1. Tells Quantum Espresso to carry out structure optimisation.   
        2. Set convergence criterion for forces in units of Ry/Bohr.   
        3. You'll need to add an `IONS` section to the input. This is the only mandatory addition.   
    !!! note "Variables of IONS"
        There are many other variables you can add for the IONS section. These include the algorithm of relaxation, adding constraints to the position of atoms, etc. If you are interested, you can look up the input description of `pw.x`.

    - Run `pw.x < CH4_opt.in > CH4_opt.out`.  
    - Let's look at the output file. The important sections are shown below. 

    ```python
        

    ```
    - Final atomic positions 
    - Final total force 
    - Final bond lengths are the same 




!!! note "Interlude"
    A short break you should take, tired if you are. 


## Ground state structure of crystals 
We have seen how the stablest structures of molecules can be found by minimising the atomic forces in DFT. In principle, it is straightforward to apply the idea for the optimisation of the structures of crystals. However, an additional consideration arises from optimising the lattice constant of the crystal alongside the atomic positions within the unit cell. This leads to the need for a few more quantities. 

### Stress in crystals
Mathematically, the stress is a rank-2 tensor, i.e. a matrix, $\sigma_{ij}$, defined through the relation  

$$
F_{i} =  \sigma_{ij}\epsilon_{j}
$$

??? tip "Tip: what does stress measure and why is it a rank-2 tensor (matrix)?"
    Stress is the ratio between the $i$-th component of the induced force (e.g. this could be the $x$-component), $F_{i}$, due to a deformation along the $j$-th direction (e.g. in the $x$- or $z$-direction), $\epsilon_{j}$. So it measures how favourable or unfavourable it is to distort the crystal along a certain direction in terms of the induced force along different directions. 

    This is also why it is a rank-2 tensor, you need to specify two pieces of information before reading the stress tensor. You need to decide which distortion direction, and which component of the induced force you would like to know from the stress tensor. 

### Pressure and bulk modulus of crystals 
Whereas stress describes the force induced on a crystal due to a distortion along a certain direction, the pressure and the bulk modulus of a crystal measures the changes in total energy due to contractions or expansions. Pressure, P, is defined as the negative of the first derivative of the PES against the volume, V. It measures the tendency of the crystal to expand or contract. Mathematically,  

$$
P = -\frac{dU}{dV}.
$$

??? tip "Tip: what does it mean to have a structure with a negative or positive pressure?"
    If the pressure is negative, the crystal structure is unstable and will contract. If the pressure is positive, the crystal structure is also unstable but it will expand. The structure is only stable when the pressure is close to zero, as first derivatives of the PES vanish at the equilibrium structure.  

The bulk modulus, $B$, per volume $V$, is defined as the negative of the second derivative of the PES against the volume. It measures the stability of the crystal against expansion or contraction. Mathematically,  

$$
B = -V\frac{d^2U}{dV^2}.
$$

??? tip "Tip: what does it mean to have a structure with a negative or positive bulk modulus?"
    If the bulk modulus is negative, the curvature of the PES against volume is positive. This means the structure is stable and it costs energy to compress or expand the crystal. If the bulk modulus is positive, then the structure is unstable and the structure will spontaneously expand or compress. The structure expands if the pressure is positive, and contracts if the pressure is negative.  

Having introduced these useful quantities, we can try to calculate these quantities using Quantum Espresso. In the following tasks, we will look at the poly(para-phenylene) (PPP) polymer. It is a chain of benzene rings which extends essentially infinitely in one direction. It is therefore a one-dimensional crystal. The unit cell of this crystal contains two benzene rings. We will first run a singe-point calculation of this crystal.  

!!! example "Task 4 - calculating stress and pressure of crystals"
    - Read the input file `PPP_strpr.in`. The important parts are shown below.  
        ```python 
         &CONTROL
            pseudo_dir = '.'
            disk_io = 'none'
            tprnfor = .true.
            tstress = .true.    #(1)!
         /
        
         &SYSTEM
            ibrav = 0
            nat = 20
            ntyp = 2
            ecutwfc = 35.0
         /
        
         &ELECTRONS
         /
        
        ATOMIC_SPECIES
         C  12.011  C.pz-vbc.UPF
         H   1.008  H.pz-vbc.UPF
        
        CELL_PARAMETERS bohr          #(2)!
        25.000000  0.000000  0.000000
         0.000000 25.000000  0.000000
         0.000000  0.000000 16.140768
        
        ...

        K_POINTS automatic        #(3)!
        1 1 3 1 1 1
        
        ```
        1. This tells Quantum Espresso to print the stress.  
        2. The CELL_PARAMETERS are specified in Bohrs.    
        3. The k-grid is specified now for a crystal calculation. We have chosen the optimal k-grid for you.   
    !!! question "Quick quiz: why do we only sample the k-points along the $z$-direction?"
        ??? success "Answer"
            Because the crystal is one-dimensional and along the $z$-direction. 
    !!! question "Quick quiz: what are the lengths of the vacuum box in the $x$ and $y$-directions?"
        ??? success "Answer"
            25 Bohrs. 

    - Run `pw.x < PPP_strpr.in > PPP.out`. 
    - Once the job has finished, look for the line saying "total   stress". You should see something like this:
    ``` python 
    
    Computing stress (Cartesian axis) and pressure
    
              total   stress  (Ry/bohr**3)                   (kbar)     P=      -25.73
      -0.00019413  -0.00001068   0.00000000          -28.56       -1.57        0.00
      -0.00001068  -0.00016867   0.00000000           -1.57      -24.81        0.00
       0.00000000   0.00000000  -0.00016197            0.00        0.00      -23.83   
    ```

    - What are the units for the stress tensor?
    ??? success "Answer"
        Ry/Bohr^3. The dimension of this unit is energy per volume, which is the same dimension as pressure. 

    - What is the total pressure? 
    ??? success "Answer"
        -25.73 kPa. Note that kPa is one of the common units for pressure.  

    - Is the lattice constant along the z-direction optimised? Why? If not, does it want to expand or contract? 
    ??? success "Answer"
        No. One evidence is that $\sigma_{zz}$ is -23.83 kPa, which is substantial. Since the sign is negative, it wants to contraact. 

    !!! note "Note"
        The stresses along $x$ and $y$ are both substantially negative. This happens when we use a large vacuum in the simulation box. But since we only have a one-dimensional crystal in here, these stress components are meaningless.   

Whereas Quantum Espresso prints the stress and pressure for you, it does not calculate the bulk modulus. You will need to distort the unit cell volume manually yourself and run multiple `pw.x` calculations to obtain the PES against the unit cell volume. Let's do this in the following task. 

!!! example "Task 5 - calculating the bulk modulus using the PES"
    - Run `python3  EXE_VOL_STRETCH.py` to generate a list of input files of different unit cell volumes. The script automatically does the generation for you so you should not need to edit anything. If you are interested, you can look at the comments in the script and see what they do 

     - Run `python3 RUN_ALL.py` to run `pw.x` on all the structures. This script is quite different from the bash scripts you have used before. If you are interested, you can read the comments too.  

    - After all the jobs are finished, run `python3  PLOT_BULK_MOD.py`. 

    - The step above should produce a graph, with a parabola fitted to it and print out a line that says "The coefficient of the parabolic fit to the PES against volume is: ..."

    - What does the PES look like? 

        ??? success "Answer"
            The graph will be uploaded soon 

    - What is the bulk modulus of this structure? 

        ??? success "Answer"
            The second derivative of a parabola with the equation $ax^2+bx+c$ is simply $a$. Therefore, the bulk modulus per volume at this structure is simply $a/V_\mathrm{struct} = .../... = $ .   

We are now on a position to carry out structural optimisation in crystals. We have seen just now that the PES plotted against unit cell volume is useful for calculating the bulk modulus. In fact, plotting the PES against different parameters of the structure is very useful for structure optimisation too.  

If you take a closer look at the PPP molecule again, you will notice a torsion angle, $\theta$, between two neighoburing rings in the unit cell. Our next task is to optimise this torsion angle using the PES.  

### Optimisation of atomic positions within the unit cell 
We are going to find the optimal torsion angle by finding the minimum in the PES, plotted against the torsion angle $\theta$. This is what we will do in the following task. 

!!! example "Task 5 - Optimising structure through PES"
    - Multiple input files named `PPP_{theta}.in` have been prepared for different torsion angles. Take a look at the values of $\theta$s used in this task. 

    - Run `bash theta_run.sh`. This script automatically runs all the `pw.x` for all the input files. It is written similar to the ones you have used for convergence tests. 

    - Once all the input files have been run, run the Python script `PLOT_PES_THETA.py`, which plots the total energy of the PPP molecule against the torsion angle.  

    - How does the PES look like and what is the optimal torsion angle? 
    ??? success "Result"
        A plot will be uploaded very soon. 

Using the PES can allow us optimise specific parameters of the structure, such as the torsion angle in the PPP chain. This is useful when we want to obtain physical insight or intuition, such as how the benzene rings interact with their nearest neighbours to arrange themselves. However, if we want to optimise the structure fully, we will have to run a "relax" calculation with Quantum Espresso. 

!!! example "Task 6 - Optimising structure through DFT relaxation"
    - Go to directory 06_PPP

    - Read the input file named `PPP_opt.in`. The torsion angle of the initial structure in this input file is 20 degrees. The input file is almost the same as the previous ones, but with the calculation type set to relax. See the snippet below: 
    ```python 
     &CONTROL
        pseudo_dir = '.'
        disk_io = 'none'
        calculation = 'relax'    #(1)!
     /
     ...  
     &IONS                  #(2)!
     /
    ```
    1. Tells Quantum Espresso to carry out relaxation. 
    2. IONS section is again required for relaxation calculation. 

    - Run `pw.x < PPP_opt.in > PPP_opt.out`. 

    - After the job finishes, open the output file in `VESTA`. As you will see, it's possible to just open the final optimized structure, or load the optimization as an animation. Since we already started quite close to the final structure, you won't see much of a change if you watch an animation. 

    - Use the `dihedral angle` option in xcrysden to find the torsion angle of the relaxed   structure. How does this compare to what you previously predicted?
    ??? success "Answer"
        The figure will be uploaded once it is ready. 


### Optimizing lattice parameters of unit cells 
Finally, we now utilise both forces and stress for optimising the lattice constant of the unit cell. We are going to relax the lattice constant of the PPP chain along the $z$-direction.  
    
!!! note "Note on choice of plane-wave cutoff"
    In doing this you should keep in mind that it can take a higher energy-cut off to converge the force or equivalently the stress as we saw at the start of this lab. Also, if you recall the number of planewaves depends on the cell volume, so if during the course of the calculation we change the volume we may also change the number of plane waves, or if we fix the number of plane waves we are changing the effective energy cut-off (the latter for Quantum Espresso, but different codes will handle this differently). So you need to be quite careful about this when optimizing lattice vectors.

!!! example "Task 7 - Optimising the unit cell"
    - Read the input file `PPP_cell.in`. You'll notice in addition to the inputs mentioned, there's also a fairly high energy cut-off, and we've lowered the SCF convergence threshold from the default. Try running this now and let's look at the output.
    ```python 
     &CONTROL
        pseudo_dir = '.'
        disk_io = 'none'
        calculation = 'vc-relax'    #(1)!
        forc_conv_thr = 1.0D-4     #(2)!
     /
     ...
     &IONS              #(3)!
     /
    
     &CELL 
     cell_dofree = 'z'    #(4)!
     /
    ``` 
        1. Tells Quantum Espresso to do a variable-cell (vc) relaxation. 
        2. Specify force tolerance. 
        3. The IONS section is needed for vc-relaxation too.  
        4. This tells Quantum Espresso to relax the $z$-component of the third lattice vector only. 

    ??? tip "Tip on constrained relaxations"
        In the input file shown just now,we have constrained the relaxation to moving one specific component of a specified lattice vector. Constraints like this is very useful especially for low-dimensional systems. If you are interested in what Quantum Espresso can can constrain for you, you can read the input description of `pw.x`.  
    - What is the final lattice constant? 

    - How can you check if the third lattice vector is fully relaxed? 

    - What is the final torsion angle? 

!!! warning "A word on the final pressure and lattice constant in the output file" 
    The final pressure at the end of the optimisation might appear large. This is because at the end of the optimization, an `scf` calculation is automatically performed starting from the optimized structure. This is needed Quantum Espresso fixes the basis set as that for the original input structure during the calculation, so if the structure has changed a lot, you may calculate a different stress when starting from the relaxed structure. You will be able to see from the final stress or pressure whether you should rerun your calculation.

    Additionally, the cell output will likely be in Quantum Espresso's representation where the cell vectors are given explicitly in Bohr along with a scaling factor `alat` which is fixed. (In our case here `alat` will be `A` converted from Angstrom to Bohr). If you want to rerun your calculation you could either input the cell using these directly, or calculate appropriate values for the input. You may need to do the latter if you want to find the new lattice constant anyway.

----------------------------------------------------------------------------------
## Summary

In this lab we have seen

- How to output forces, stresses, and pressures in Quantum Espresso calculations.  

- How to check if the forces have converged.

- How the PES allows to estimate the bulk modulus and optimise specific structural parameters.

- How to use these forces in a calculation to optimize the atomic positions, where they are moved until the forces on the atoms are less then some threshold value.

- How the stresses can be used to optimize the structure, where the lattice constant is changed until the stress on the system is less than some threshold value.

