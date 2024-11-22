Structural Optimisation
===============================
For this week and the next, we will focus on predicting the structural properties of materials. These include finding the stablest structure of molecules and crystals, and calculating the vibrational patterns in molecules and crystals.

In this week, we will learn how to calculate the key quantities for describing the stability of the structures of molecules and crystals. You will have learnt them from the lectures, and we will give a brief review for them in this lab before the calculations to refresh your memory. We will begin by briefly reviewing the very important potential energy surface (PES). 

## The potential energy surface
The potential energy surface (PES) is the surface obtained by plotting the total energy of a molecule or crystal against different atomic positions. It is commonly denoted as $U(\{\mathbf{R}\})$, where $\mathbf{R}$ represents the set of atomic positions of all atoms. Many useful quantities associated with a given structure, such as the atomic forces in a molecule, the stress and pressure in a crystal unit cell can be calculated by evaluating the relevant derivatives of the PES.  

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

In following tasks, we will go through how the atomic forces are calculated using Quantum Espresso. We have prepared a distorted methane molecule, $\mathrm{CH_4}$, where one of the hydrogen atoms (the one sitting on the z-axis above the $\mathrm{C}$ atom) is pushed closer to the carbon atom than others. 

!!! example "Task 1 - Examining atomic forces in the output file"
    - Go to the directory `01_ch4_scf`. 
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
    - What is the size of the supercell used in this calculation? 
    ??? success "Answer"
        20 Angstrøms. 

    - Why are these atomic forces expected for this distorted methane molecule? 
    ??? success "Answer"
        Since the $\mathrm{H}$ atom (atom 2) above the $\mathrm{C}$ atom (atom 1) is approaching the $\mathrm{C}$ atom, the covalent $\mathrm{C}-\mathrm{H}$ bond between these two atoms is compressed. The two atoms will repel each other to uncompress the bond. Atom 2 should experience an upward force and atom 1 should experience a downward force. This is indeed the case in the output file.  


    - What are the units of atomic forces in Quantum Espresso? 
    ??? success "Answer"
        Ry/Bohr. (1 Ry/Bohr is equivalent to 25.7 eV/Å, check this yourself!). In this task, the $z$-component of the force acting on the $\mathrm{H}$ atom is 0.0527 Ry/Bohr. In the linear approximation, this means drifting the $\mathrm{H}$ atom in the positive $z$-direction by 0.1 Bohr should lead to a reduction in total energy of about 5.27 mRy. This is equal to 0.071 eV, which is substantial. This is because the covalent bond is quite compressed. 

!!! Note 
    The `Total force` listed here is the square root of the sum of _all_ of the force components squared rather than the sum of the magnitudes of the individual forces on the atoms. It is likely because the number is intended more as a guide to check overall accuracy. If the `Total SCF correction` is comparable to the `Total force` it usually means you need to try to better converge the SCF cycle (via `conv_thr`).

### Convergence test 
Just as we have to check the convergence of the total energy against the PW cutoff, we have to check the convergence of atomic forces against the PW cutoff too. Let's do this in the next task.  

!!! example "Task 2 - Examining convergence"
    <!--In this task, multiple input files of the same distorted methane molecule are prepared. Each input file is labelled by their PW cutoff.--> 
    
    - Go to the directory `02_ch4_convergence`. 

    - Similar to what you have done in Lab 3, make 9 copies of the input file `CH4_20.in` and name them `CH4_i.in`, where `i` runs from 30 to 100 in steps of 10. Edit each input file so that the value of `ecutwfc` is `i` Ry. 

    <!--Run `python3 ecutwfc_file_build.py` script to produce a list of input files. This script is very similar to the one used in Lab 3. This will create input files with PW cutoffs from 10 to 100 Ry. --> 
    
    - Run `pw.x` to perform a DFT calculation for each input files.
    - Create a file named `data.txt` that records the values of PW cutoff in the first column, and the values of `Total force` in the second column. 
    <!-- - Run `bash ecutwfc_force.sh` to run all the `pw.x` calculations as in Lab 3.--> 

    - Which quantities are being plotted by the script for determining the convergence? 

    ??? success "Answer" 
        The script plots the fractional difference of the total force with respect to the best-converged total force $(F^{\mathrm{tot}}_{\mathrm{best}}-F^{\mathrm{tot}}_{i})/F^{\mathrm{tot}}_{\mathrm{best}}$.   
    
    - How does the convergence curves look like? 
    ??? success "Answer"   

    !!! note "Note: rate of convergence of forces"
        Atomic forces converge slower with respect to PW cutoff compared to the total energy. So if your goal is to find the stablest structure, make sure the total force is converged in addition to the total energy.  

### Optimisation of molecular geometry 
Now that we have seen how to calculate the atomic forces associated with a particular structure, we can use Quantum Espresso to optimise the structure of molecules.

!!! example "Task 3 - Optimising the molecular structure"
    - Go to the directory `03_ch4_opt`. 
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

    - Run `pw.x < CH4_opt.in > CH4_opt.out`.  
    - In the output file, navigate to the line `End of BFGS Geometry Optimization`. The important sections above the end of the output file are shown below. 
    ```python
     ...

     number of scf cycles    =   9 #(1)!
     number of bfgs steps    =   8 #(2)!
     
     ... 

     Forces acting on atoms (cartesian axes, Ry/au):

     atom    1 type  1   force =    -0.00000070    0.00000000   -0.00008813   #(3)! 
     atom    2 type  2   force =     0.00000198    0.00000000   -0.00008060
     atom    3 type  2   force =    -0.00003361    0.00000000    0.00005439
     atom    4 type  2   force =     0.00001617    0.00002899    0.00005717
     atom    5 type  2   force =     0.00001617   -0.00002899    0.00005717

     Total force =     0.000165     Total SCF correction =     0.000306 #(4)! 
       
    ... 

     bfgs converged in  10 scf cycles and   8 bfgs steps #(5)! 
        
    ```
        1. The latest number of scf cycles is shown after each scf cycle. 
        2. The latest number of relaxation steps is shown after each relaxation step. Here, a relaxation step is called a "bfgs step" in our calculation. This is just the name of our relaxation algorithm.
        3. The forces on each atom is almost zero. 
        4. Note that the total forces are below `forc_conv_thr`.
        5. The total number of scf cycles is 10, and the number of relaxations steps taken is 8. Note that there tends to be more scf cycles than relaxations steps in geometrical optimisation.  
    

    ??? success "What is the equilibrium bond length?"
        1.08 Å. The experimental value is 1.09 Å so this is in good agreement. You can check that all the bond lengths are the same.  


    !!! note "Variables of IONS"
        There are many other variables you can add for the IONS section. These include the algorithm of relaxation, adding constraints to the position of atoms, etc. If you are interested, you can look up the input description of `pw.x`.



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

??? tip "Tip: what does it mean to have a structure with a negative bulk modulus?"
    If the bulk modulus is negative, the curvature of the PES against volume is positive. This means the structure is close to being stable. The structure expands if the pressure is positive, and contracts if the pressure is negative.  

Having introduced these useful quantities, we can try to calculate these quantities using Quantum Espresso. In the following tasks, we will look at carbon diamond. The unit cell of this crystal contains two carbon atoms. We will first run a single-point calculation of this crystal.  

!!! example "Task 4 - calculating stress and pressure of crystals"
    - Read the input file `CD_scf.in`. The important parts are shown below.  
        ```python 
         &CONTROL
            pseudo_dir = '.'
            prefix = 'CD'
            disk_io = 'low'
            tprnfor = .true.
            tstress = .true.    #(1)!
         /
        
        ... 

        ATOMIC_POSITIONS crystal   #(2)!
         C 0.00 0.00 0.00
         C 0.25 0.25 0.25
        
        K_POINTS automatic        #(3)!
        8 8 8 0 0 0
        
        ```
        1. This tells Quantum Espresso to print the stress.  
        2. The fractional coordinates are used for atomic positions. This allows Quantum Espresso to speed up calculations with symmetry detection.     
        3. The k-grid is specified now for a crystal calculation. We have chosen the optimal k-grid for you.   

    - Run `pw.x < CD_scf.in > CD.out`. 
    - Once the job has finished, look for the line saying "total   stress". You should see something like this:
        ```python 
             Computing stress (Cartesian axis) and pressure
        
                  total   stress  (Ry/bohr**3)                   (kbar)     P=      -14.69
          -0.00009989  -0.00000000  -0.00000000          -14.69       -0.00       -0.00
          -0.00000000  -0.00009989  -0.00000000           -0.00      -14.69       -0.00
           0.00000000   0.00000000  -0.00009989            0.00        0.00      -14.69
        ```
    
    - What are the units for the stress tensor?
    ??? success "Answer"
        Ry/Bohr^3. The dimension of this unit is energy per volume, which is the same dimension as pressure. 

    - What is the total pressure? 
    ??? success "Answer"
        -14.69 kPa. Note that kPa is one of the common units for pressure.  

    - Is the lattice constant optimised? Why? If not, does it want to expand or contract? 
    ??? success "Answer"
        No. One evidence is that $\sigma_{zz}$ is -14.69 kPa, which is substantial. Since the sign is negative, it wants to contract. 

    !!! note "Note"
        Sometimes, when simulating molecules (or 2D crystals) with a large supercell, the stresses can substantially negative. This happens when we use a large vacuum in the simulation box. So the components of the stress along the directions of vacuum is large. These stress components are meaningless.   

Whereas Quantum Espresso prints the stress and pressure for you, it does not calculate the bulk modulus. You will need to distort the unit cell volume manually yourself and run multiple `pw.x` calculations to obtain the PES against the unit cell volume. Let's do this in the following task. 

!!! example "Task 5 - calculating the bulk modulus using the PES"
     - Go to the directory `05_cd_bm`. You should find an input file called `CD_com.in`. This is the reference input file representing a diamond structure with an unrelaxed lattice parameter $A_0=3.502Å$. We are going to change the volume of the unit cell by changing the lattice parameter. In other words, we are going to compress or elongate the lattice vectors to produce the desired change in volume.  
     
    - Make 10 copies of this input files named `CD_i.in`, where `i` run from -0.08 to 0.08 in steps of 0.02. Here `i` denotes the fractional change in volume, $\Delta V/V_0$, where $V_0$ is the original volume, and $\Delta V$ is the change in volume. In other words, if the new volume is $V'$, then $\Delta V = V'-V_0$. Note that one of the files should be `CD_0.00.in`, which is the input file having the same cell volume as the reference input file. 

    - Read the input file `CD_com.in`. Notice how we have to carry out a relaxation calculation for each cell volume to allow the atoms to move to the correct positions for each volume. 

    - Now, for each fractional change in volume, calculate the required strain $\varepsilon$ that needs to be applied on the lattice parameter to produce such a change in volume. The correct equation for calculating the strain is derived in the following box. 

    ??? success "What is the equation relating the strain $\varepsilon$ and fractional change in volume $\Delta V/V$?"
        One can represent the fractional change in volume as $\alpha$. So that the volume changes as $V_0 \rightarrow V'=(1+\alpha)V_0$, where $V_0$ is the original volume. 

        The strain is defined as the fraction $\varepsilon$ such that the change in lattice parameter is $A_0 \rightarrow A'=(1+\varepsilon)A_0$, where $A_0$ is the original lattice parameter. 

        Since $V$ is always proportional to $A^3$, we have
        $$
        V'/V_0 = (A'/A_0)^3,
        $$
        so that
        $$
        (1+\alpha) = (1+\varepsilon)^3.
        $$
        Therefore, 
        $$
        \varepsilon = (1+\alpha)^{1/3}-1. 
        $$
        Hence, the new lattice parameter is given by 
        $$
        A' = (1+\alpha)^{1/3}\times A_0
        $$
        
        So given a desired fractional change in volume $\alpha_i$, you should calculate the corresponding new lattice parameter $A_i'$ using the last equation. 
    
    - Edit the value of the lattice parameter `A` for each input file using the values you have calculated based on the last equation. We recommend using Excel or Python to calculate these values.  
    
    <!-- - Run `python3  EXE_VOL_STRETCH.py` to generate a list of input files of different unit cell volumes. The script automatically does the generation for you so you should not need to edit anything. If you are interested, you can look at the comments in the script and see what they do -->
    <!-- - Run `python3 RUN_ALL.py` to run `pw.x` on all the structures. This script is quite different from the bash scripts you have used before. If you are interested, you can read the comments too. -->

    - Run `pw.x` on each input file, and create a file named `data.txt`, with the first column containing all the cell volumes $V$, and the second column containing the total energies of the corresponding scf calculation. 

    - Run `python3 BULK_MOD.py`. 

    - The step above should produce a graph, with a parabola fitted to it and print out a line that says "The coefficient of the parabolic fit to the PES against volume is: ..."

    - What does the PES look like? 

        ??? success "Answer"
            The graph will be uploaded soon 

    - How is the bulk modulus calculated? What is its value for this structure? 

        ??? success "Answer"
            The second derivative of a parabola with the equation $ax^2+bx+c$ is simply $2a$. Therefore, the bulk modulus per volume at this structure is simply $-2aV_{0} = ... $ .   

    - Is the structure stable? If not, does it want to expand or contract? 

        ??? success "Answer"
            The bulk modulus is negative, so the structure is close to being stable. Since the pressure is positive, the unit cell wants to expand. In other words, the reference lattice constant is smaller than the equilibrium lattice constant. 

    - Read the output file for `i`=0.08. How many relaxation steps does it take for the calculation to finish? 

        ??? success "Answer"
            None. In this particular setup, through the use of fractional coordinates, the carbon atoms are guaranteed to sit at the high-symmetry positions, so that the net force acting on them vanishes. So for any given volume, the atomic positions are already relaxed. In fact, if you read any other output file, the number of bfgs steps should also be 0.  
            

We will now carry out the structural optimisation in carbon diamond to find the equilibrium lattice parameter.  

### Optimizing lattice parameters of unit cells 
Finally, we now utilise both forces and stress for optimising the lattice constant of the unit cell. We are going to relax the lattice constant of carbon diamond.  
    
!!! note "Note on choice of plane-wave cutoff"
    In doing this you should keep in mind that it can take a higher energy-cut off to converge the force or equivalently the stress as we saw at the start of this lab. Also, if you recall the number of planewaves depends on the cell volume, so if during the course of the calculation we change the volume we may also change the number of plane waves, or if we fix the number of plane waves we are changing the effective energy cut-off (the latter for Quantum Espresso, but different codes will handle this differently). So you need to be quite careful about this when optimizing lattice vectors.

!!! example "Task 6 - Optimising the unit cell"
    - Go to the directory `06_cd_opt` and read the input file `CD_opt.in`. You'll notice in addition to the inputs mentioned, there's also a fairly high energy cut-off, and we've lowered the SCF convergence threshold from the default. Run this now with `pw.x`.
    ```python 
     &CONTROL
        pseudo_dir = '.'
        prefix = 'CD'
        disk_io = 'low'
        calculation = 'vc-relax'    #(1)!
        etot_conv_thr = 1.0D-6
        forc_conv_thr = 1.0D-4     #(2)!
     /
     ...
     &IONS              #(3)!
     /
    
     &CELL     #(4)!
     /
    ``` 
        1. Tells Quantum Espresso to do a variable-cell (vc) relaxation. 
        2. Specify force tolerance. 
        3. The IONS section is needed for vc-relaxation too.  
        4. The CELL section is also needed for vc-relaxation.  

    Now inspect the output file. It is very similar to the one for relaxing molecules in task 3.  

    - What is the final lattice constant? 
    ??? success "Answer" 

    - How can you check if the lattice parameter is fully relaxed? 
    ??? success "Answer" 
          The pressure is very small. 

 <!--    ??? tip "Tip on constrained relaxations"
        In the input file shown just now,we have constrained the relaxation to moving one specific component of a specified lattice vector. Constraints like this is very useful especially for low-dimensional systems. If you are interested in what Quantum Espresso can can constrain for you, you can read the input description of `pw.x`. -->  

!!! warning "A word on the final pressure and lattice constant in the output file" 
    The final pressure at the end of the optimisation might appear large. This is because at the end of the optimization, an `scf` calculation is automatically performed starting from the optimized structure. This is needed Quantum Espresso fixes the basis set as that for the original input structure during the calculation, so if the structure has changed a lot, you may calculate a different stress when starting from the relaxed structure. You will be able to see from the final stress or pressure whether you should rerun your calculation.

    Additionally, the cell output will likely be in Quantum Espresso's representation where the cell vectors are given explicitly in Bohr along with a scaling factor `alat` which is fixed. (In our case here `alat` will be `A` converted from Angstrom to Bohr). If you want to rerun your calculation you could either input the cell using these directly, or calculate appropriate values for the input. You may need to do the latter if you want to find the new lattice constant anyway.

----------------------------------------------------------------------------------
## Summary

In this lab we have seen:

- How to calculate atomic forces, stress, and pressure in Quantum Espresso. 

- How to check if the forces have converged.

- How the PES allows us to estimate the bulk modulus and pressure.

- How to optimise the structure of molecules, atomic positions in a crystal, and the lattice parameter of a crystal. 
