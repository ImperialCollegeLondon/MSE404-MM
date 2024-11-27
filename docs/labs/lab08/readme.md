# Vibrational Normal Modes and Phonons
This week we continue studying the structural properties of molecules and crystals. We will calculate the vibrational properties of molecules and solids. 

As you will have seen in the lectures, the vibration patterns and the vibrational energies can also obtained from the potential energy surface (PES). More specifically, the PES allows us to derive an equation of motion of the vibrations for the atoms around the equilibrium positions in the equilibrium structure. Importantly, in deriving the equation of motion, two important matrices arise: the force constant matrix and the dynamical matrix. 

Similar to last week, we will give a brief review of these matrices and when they become useful, leaving the mathematical details to the lectures. Then we will learn how to calculate them using Quantum Espresso. Let's now begin our dive into vibrations in molecules. 

## Vibrations in molecules: basic theory
When molecules are excited (e.g. by shining light or raising temperature), they gain energy and the atoms starting vibrating around their equilibrium position. The atoms vibrate collectively, raising the energy of the molecule. The whole vibrational analaysis boils down to analysing this change in total energy.  

### Force constant matrix 
The change in total energy around the equilibrium structure due to small displacements, $\Delta U$, can be written as,  
$$
\Delta U = \sum_{I\alpha J\beta} K_{I\alpha J\beta}u_{I\alpha}u_{J\beta},
$$
where $u_{I\alpha}$ denote the $\alpha$-component of the displacement of the atom labelled by $I$, and similarly for $u_{J\beta}$. The matrix $K_{I\alpha J\beta}$ is called the force constant matrix. It is defined by the equation 
$$
K_{I\alpha J\beta} = \frac{\partial^2 U}{\partial u_{I\alpha}\partial{u_{J\beta}}}.
$$
Note that $I$ and $J$ both run from $1$ to $N$, where $N$ is the number of atoms, and $\alpha$ and $\beta$ run from $1$ to $3$, with the convention $x\rightarrow1, y\rightarrow2, z\rightarrow3$.

??? tip "Tip on understanding the indices of $\mathbf{K}$"
    We have to understand the indexing of the matrix $\mathbf{K}$ carefully by looking at its physical meaning. It is called the force constant matrix because the entries physically mean the $\alpha$-component of the induced atomic force acting on atom $I$ due to the atom $J$ moving along direction $\beta$ for one unit. This is also why it has four indices: four pieces of information have to be specified before calculating the force.

??? tip "Tip on writing down $\mathbf{K}$ as a two-dimensional matrix"
    The matrix $\mathbf{K}$ can be written as a two-dimensional matrix although its entries are labelled by four indices. We can combine the indices into two pairs: $\{I\alpha\}$ and $\{J\beta\}$. In this way we can write $\mathbf{K}$ as a $3N\times3N$ matrix, where $N$ is the number of atoms and the factor of 3 comes from the fact that there are three Cartesian directions, $\{x,y,z\}$.

??? question "Optional quick quiz - for your own understanding"

    !!! question "What does the entry $K_{1123}$ physically mean?" 
        ??? success "Answer"
            It is the force acting on atom 1 along the $x$-direction when atom 2 (and only atom 2) moves along the $z$-direction 
    
    !!! question "What does the entry $K_{1111}$ physically mean?"
        ??? success "Answer"
            It is the $x$-component of the force experienced by atom 1 when it moves along the $x$-direction. In most cases this should be positive because atom 1 should experience a restoring force that drives it back to the equilibrium position.  
    
    !!! question "What geometrical features of the PES do the entries $K_{1111}$ and $K_{1123}$ mean?"
        ??? success "Answer"
            $K_{1111}$ represents the curvature of the PES at the equilibrium structure along $x_1$, where $x_1$ is the $x$-coordinate of atom 1. Similary, $K_{1123}$ represents the curvature at the equilibrium structure when both $x_1$ and $z_2$, where $z_2$ is the $z$-coordinate of atom 2, are varied around the equilibrium. 

<!--!!! note "The force constant matrix is a symmetric matrix"-->

The force constant matrix is useful for calculating the total energy change due to any displacement pattern. However, it does not allow us to calculate the dynamical properties, specifically the vibration patterns and their frequencies. 

### Dynamical matrix 
To calculate the dynamical properties, we need one more matrix. To get there we have to combine the above equations with Newton's second law, as shown in the lectures. By doing this another useful matrix appears in the equation of motion, and that matrix is therefore called the dynamical matrix.  

The equation of motion for the vibrations, which has been derived for you in the lectures, is,
$$
\omega^2\mathbf{v}=\mathbf{D}\mathbf{v}, 
$$
where $v_{I\alpha}=\frac{1}{\sqrt{M_I}}u_{I\alpha}$ are the mass-weighted displacements, $\omega$ denotes the vibration frequencies, and the matrix $\mathbf{D}$ is the dynamical matrix. It is straightforward to show that the entries of $\mathbf{D}$ are related to those of $\mathbf{K}$ by 

$$
D_{I\alpha J\beta}=\frac{1}{\sqrt{M_I M_J}}K_{I\alpha J\beta}. 
$$

??? tip "Tip: what does the entries of $\mathbf{D}$ physically mean?"
    The physical meaning of the entries of $\mathbf{D}$ could be made more transparent by comparing the equation of motion to that of a one-dimensional classical spring-mass oscillator, with spring constant $k$ and mass $m$, which is 
    $$
    \omega^2x=\frac{k}{m}x,
    $$
    where $x$ is the displacement and $\omega$ is the oscillation frequency. By comparing this equation to $\omega^2\mathbf{v}=\mathbf{D}\mathbf{v}$, we see that $\mathbf{D}$ plays the same role as $k/m$ - it measures the ratio of the strength of the spring constant of the atomic forces to the mass of the vibrating object. 

    In other words, it measures how difficult it is to move the atoms of masses $M_I$ due to the atomic forces induced by the displacement. In fact, the units of the entries of $\mathbf{D}$ have the same dimension as that of $k/m$. 

??? question "Optional quick quiz - for your own understanding"

    !!! question "Which pair of atomic displacements are being considered for $D_{1x2y}$?" 
        ??? success "Answer"
             The $x$-displacement of atom 1 and the $y$-displacement of atom 2. 

    !!! question "Which pair of atomic displacements are being considered for $D_{3z1y}$?" 
        ??? success "Answer"
             The $z$-displacement of atom 3 and the $y$-displacement of atom 1.

<!--!!! note "The dynamical matrix is also a symmetric matrix."-->

The vibration frequencies of molecules are quantised. The vibration patterns that corresponds to these frequencies are called the __normal modes__ of vibrations. All that remains is to find the vibration patterns and frequencies of the molecule.  

### Vibrational normal modes
The normal modes of vibrations are the eigenvectors of the dynamical matrix $\mathbf{D}$, as discussed in the lectures. The corresponding eigenvalues are their vibration frequencies. These normal modes tend to be symmetrical. In other words, the displacement patterns tend to demonstrate all or some of the symmetries of the molecules, as we shall see below.  

The number of normal modes of any molecule can be calculated easily. If there are $N$ atoms in the molecule, there will be $3N$ normal modes. 

## Vibrations in molecules: Quantum Espresso (overview) 
Density functional theory can be used for calculating the dynamical matrix $\mathbf{D}$ of molecules. This requires the use of density functional perturbation theory (DFPT), which you have been introduced to in the lectures. Now we can try to run our first normal mode calculation with Quantum Espresso. 

!!! warning "Heads-up"
     The way Quantum Espresso finds normals modes of molecules invokes extra rules and theories. These theories are beyond the scope of this course, but their use in Quantum Espresso means that there are a number of items you have to be careful of when preparing input files and understanding the output files. 

    In this lab, we will cover the necessary items which you **must** include in the input files, otherwise the calculation will not give you physical results. We will only outline briefly why these items have to be added, and leave the theoretical details as extra notes for the interested (and the mathematically bold!).  

## Vibrations in molecules: Quantum Espresso (calculation)
We will walk though how to find the normal modes of a methane molecule using Quantum Espresso in the following tasks. You will learn how to use an additional module called `ph.x`, for calculating the normal modes of a molecule. This is a two-step calculation. 

### Step 1: run the `pw.x` calculation

The first step is to carry out a single-point calculation using the equilibrium structure of methane. 
!!! example "Task 1a - run `pw.x`"
    - Read the input file `01_CH4_scf.in`. The variables are still the usual ones that you have seen in the previous labs. 

    - Note however that we have explicity ask for an unshifted 1x1x1 k-grid, which simply means we are asking for a $\Gamma$-point calculation. This is written to allow Quantum Espresso to use some optimisation routine. 
    
    - The position of the $\mathrm{H}$ atoms have optimised and rotated so that the magnitude of their $x,y,z$-coordinates are the same. Again this is for allowing Quantum Espresso to pull some optimisation tricks. 
    
    - Now run `pw.x < 01_CH4_scf.in > pw.out`.

### Step 2: run the `ph.x` calculation
The second step is to carry out the DFPT calculation using `ph.x`. Two additional items appear when you go through the input file for `ph.x`. We will explain them as you go.  

!!! example "Task 1b - run `ph.x`"
    - A second input file for calling the additional `ph.x` module is needed. This has been prepared and is named `02_CH4_ph.in`. Let's take a look at it: 
    ```bash 
        phonons of CH4 (gamma only)  #(1)!
         &INPUTPH   #(2)!
          tr2_ph = 1.0d-15  #(3)!
          asr = .true.   #(4)!
         /
        0.0 0.0 0.0    #(5)!
    ```
        1. The first line can be any informative comment that helps yourself. 
        2. This is a new card required for performing `ph.x` calculations
        3. This is the scf convergence criterion for this normal mode calculation. Notice how it is much stricter than usual.  
        4. This is the first new item. 'asr' here stands for the Acoustic Sum Rule. The acoustic sum rule prevents the frequencies of the normal modes from going negative, which most of the time are erroneous numerical artefacts. This will be explained at the end of the lab if you are interested. 
        5. This is the second new item. It restricts the normal mode calculation to the $\Gamma$-point. We will explain what this means later in the lab. But you are encouraged make a guess!  

    - If you are interested, refer to the `INPUT_PH.txt` file of the quantum espresso documentation. Alternatively, you can look up   [Phonon](https://www.quantum-espresso.org/Doc/INPUT_PH.html)

    - Run `ph.x < 02_CH4_ph.in > ph.out`. 

    - You have successfully finished your first normal mode calculation!  

### Output files of ph.x 
There are two major output files from a `ph.x` calculation. The first one is the `.out` file, which will tell you information about the DFPT calculations. It is not the most important for understanding the normal modes of the methane molecule so we will leave it out. The second file is the `matdyn` file which contains the dynamical matrix. This is the most important file for finding the normal modes.  

!!! example "Task 2 - Reading the dynamical matrix from Quantum Espresso"
    - The most important part for us begins with ` Dynamical  Matrix in cartesian axes`. The important sections are shown below.  


    - Let's go through this file bit by bit.  
    ```python
         q = (    0.000000000   0.000000000   0.000000000 ) #(1)!
    
         ... 

        1    2  #(2)! 
     -0.24917799   0.00000000    -0.17192431   0.00000000    -0.17192431   0.00000000  #(3)! 
     -0.17192431   0.00000000    -0.24917799   0.00000000    -0.17192431   0.00000000  #(4)!
     -0.17192431   0.00000000    -0.17192431   0.00000000    -0.24917799   0.00000000  #(5)!
         
         ... 
    
         Diagonalizing the dynamical matrix
    
         q = (    0.000000000   0.000000000   0.000000000 ) 
    
     ************************************************************************** #(6)! 
         freq (    4) =       2.365513 [THz] =      78.905014 [cm-1] #(7)!
     (  0.000000  0.000000  0.000092  0.000000  0.000000  0.000000 ) #(8)!
     ( -0.353554  0.000000  0.000093  0.000000  0.353553  0.000000 ) #(9)!
     ( -0.353554  0.000000  0.000093  0.000000 -0.353553  0.000000 ) 
     (  0.353554  0.000000  0.000093  0.000000  0.353553  0.000000 ) 
     (  0.353554  0.000000  0.000093  0.000000 -0.353553  0.000000 ) 
    ```
        1. Tells us that this is a vibrational $\Gamma$ point calculation. 
        2. The matrix shown below is going to be $D_{1{\alpha}2{\beta}}=\frac{1}{\sqrt{M_{1}M_{2}}}K_{1{\alpha}2{\beta}}$. I.e. the dynamical matrix associated with the displacements of atoms 1 and 2 along Cartesian directions $\alpha$ and $\beta$. 
        3. These entries are $D_{1x2x}$, $D_{1x2y}$, and $D_{1x2z}$ respectively. Note that the entries on the second, fourth, and sixth columns are the imaginary parts of the entry.  
        4. These entries are $D_{1y2x}$, $D_{1y2y}$, and $D_{1y2z}$ respectively. 
        5. These entries are $D_{1z2x}$, $D_{1z2y}$, and $D_{1z2z}$ respectively. 
        6. The entries of the diagonalised dynamical matrix will be found below this line. 
        7. The vibration frequency of the seventh normal mode. 
        8. The displacement of each atom in the seventh normal mode. This line says that the displacement of the C atom (i.e. atom 1) is (0,-0.447,0). Note that the numbers on the second, fourth, and sixth columns are the imaginary parts of the displacements. So if these numbers are not zero, you need to check your calculation. 
        9. Similar to the previous line, this line says that the displacement of the H atom (i.e. atom 2) is (0,-0.447,0). 

    - So everything you want to know about the normal modes of a molecule will be contained in the `matdyn` file. Answer the following questions by reading this file. 

    ??? success "How many normal modes are there? Is this what we should expect?"
        15, since there are 15 frequencies in the file. This is expected because there are 5 atoms so the number of normal modes is $5\times3=15$.  

    ??? success "What are the **distinct** frequencies (up to 2 decimal places)?"
        The distinct frequencies are 0.01 THz, 2.37 THz, 37.37 THz, 44.29 THz, 87.07 THz, and 90.89 THz. Note that there are always three normal modes with nearly zero frequencies. 

    ??? success "What are the degeneracies of each distinct frequency?"
        In the order of increasing energy, it is 3, 3, 3, 2, 1, 3. It is typical to find degenerate energy levels in symmetrical molecules.  

https://people.chem.ucsb.edu/laverman/GauchoSpace/methane_vib.html

??? note "Notes - why do the three lowest normal modes always have zero frequency?"
    These normal modes corresponds to the rigid translations of the molecule, i.e. drifting the molecule without distorting any chemical bonds. Therefore, there should not be any molecular vibrations, hence the zero vibrational frequency. You can check that for the first three normal modes, the displacements of all atoms are the same, implying rigid translation of the molecule. In fact, the three modes are simply rigid translations along the $y$, $x$, and $z$-directions respectively (check this yourself!).   

<!--??? notes "Extra notes: reading `ph.out` file" -->
<!--    The `ph.out` file contains quite a few information. Most of them are for your reference only and will be more useful for vibrational spectroscopy experiments. The following is what you will find: -->
    
<!--    1.  The frequencies of each normal mode in THz and $\mathrm{cm^{-1}}$. This is the most important for us. This will be written for example as `freq (    7) =      37.362598 [THz] =    1246.282111 [cm-1]`. -->
<!--    2.  The iterations performed for the DFPT calculations. This will be useful for monitoring and improving convergence of the calculations. -->
<!--    3.  The symmetry of each normal mode. This is more relevant for vibrational spectroscopy. This will be written for example as `freq (   7-   9) =       1246.3  [cm-1]   --><!-- T_2  G_15 P_4   I+R`. -->
    
<!--        ??? tip "Extra notes: what does the symbols `T_2  G_15 P_4   I+R` mean? (For your interest only)" -->
<!--            The first symbol tell us the degeneracy of this normal mode. The symbols `G_15 P_4` is the label of this normal mode in a language familiar to vibrational spectroscopists. Lastly, the symbol `I+R` tells us that both **I**nfrared and **R**aman spectroscopy can detect this normal mode.  -->
??? tip "Extra notes - symmetry reduction in Quantum Espresso (for the interested)"
     The purpose of symmetry reduction in Quantum Espresso is to minimise the number of calculations required. In the case of molecular normal mode calculations, it is about minimising the number of DFPT calculations. This is a powerful procedure because many molecules demonstrates symmetry. The code is much better able to detect the symmetry of the molecule. This is very important for calculations of the vibrations. If the code understands that all the hydrogen atoms are equivalent by symmetry, it only needs to calculate the derivatives of the energy with respect to one of the hydrogen positions, and then it can populate the full dynamical matrix based on the symmetry of the system.
     Quantum Espresso tells you what symmetries has been identified in the molecule in `ph.out` in a rather specific notation. It is written in a convention more familiar to the vibrational spectroscopy community. You can look up Raman spectroscopy and Mulliken symbols if you are interested.
        
    
<!--    4.  The number of normal modes identified, and some information on how the frequency of that normal mode is calculated. This is really for reference only. You will see something like `Representation     7      1 modes - Calculated using symmetry` -->

### Another example: the CFC molecule
Now that you have run your first normal mode calculation, let's practice it on another molecule called CFC. It is a major green gas molecule. It has a similar structure as the methane molecule, but the H atoms are replaced by three Cl and one F atoms. 

!!! example "Task 3 - normal modes of the CFC molecule"
    - Go to the directory ... . There should be two input files. One for `pw.x` and one for `ph.x`. Note that in our input file, the first atom is the C atom, the second to fourth atoms are the Cl atoms, and the fifth atom is the F atom. 

    - Run `pw.x` on ... and `ph.x` on ... 

    ??? success "What is the highest frequency of the normal modes?"
        The highest frequency mode is ... , 

    ??? success "Which three modes are dominated by the vibration of the CF bond?"
        Modes 13, 14, and 15. In these modes, the displacements of the F atom are much larger than the the Cl atoms. So these are the modes where the vibration of the CF bond dominates the normal mode. 

    ??? success "Is the highest vibrational frequency of the CFC molecule lower or higher than that of the methane molecule? Why is this the case?" 
        It is lower than that of the methane molecule. It corresponds to the normal mode where the C-F bond dominates the vibration. We can do a back-of-the-envelope estimation of this oscillation frequency. The oscillation frequency of the C-F bond, $\omega_{\mathrm{CF}}$, can be estimated from the oscillation frequency of the C-H bond, $\omega_{\mathrm{CH}}$, with the equation $\frac{\omega_{\mathrm{CF}}}{\omega_{\mathrm{CH}}}=\sqrt{\frac{k_\mathrm{CF}}{k_\mathrm{CH}}\frac{m_H}{m_F}}$ . The $k$'s denotes the spring constant associated with each chemical bond.
        
        We can replace the ratio of spring constants with the ratio of bonding energies of each bond for a crude estimate. Here $E_{\mathrm{CH}}=$, and $E_{\mathrm{CF}}=$ (from Wikipedia), and the atomic masses of H and F atoms are 1 and 19 respectively. Using $\omega_{\mathrm{CH}}=$ from the previous task we obtain $\omega_{\mathrm{CF}}=$ . This is close to the frequencies of the 13th, 14th and 15th modes.  


## Vibrations in crystals: basic theory
We have seen how vibration patterns in molecules are quantised in energy and possess similar symmetries to that of the molecule. What about vibrations in crystals?

Vibrations in crystals are simliar to those in molecules with one important difference. First of all, the atoms in each unit cell essentially behaves like a "molecule". They therefore vibrate along specific normal modes too. But the one important difference between crystals and molecule is that the vibrations between neighbouring unit cells **couple**. This leads to a collective vibration where atoms from multiple unit cells vibrate along the same normal mode, but the amplitude of vibration from one unit cell to another is modulated by a wave packet. Essentially we get a "Mexican wave" of vibrations across multiple unit cells. Such a collective vibration of atoms across multiple unit cells is called a __phonon__. 

### Phonons 
A phonon is a collective wave of atomic vibrations over multiple unit cells. Each phonon is described by two labels: the branch $\nu$ and its wavevector $\mathbf{q}$. The branch refers to the unit cell normal mode along which all atoms vibrate. The wavevector describes the propagation direction of the phonon and its spatial periodicity. The magnitude of the wavevector is given by
$$
q=\frac{2\pi}{\lambda},
$$
where $\lambda$ is the wavelength of the phonon. 

Similar to electrons (or any other quantum mehcanical particle), the momentum of the phonon is simply $\hbar\mathbf{q}$. So some textbooks will also call $\mathbf{q}$ the momentum of the phonon for simplicity.  

!!! warning "Using different symbols for the wavevectors of electrons and phonons"
    Starting from now, we will use $\mathbf{q}$ to denote the wavevector of a phonon and $\mathbf{k}$ to denote the wavevector of an electron to avoid confusion. 

### Phonon band structure 
While the branch and wavevector of a phonon tell us about the spatial pattern of the atomic vibrations, they do not tell us about the temporal features, namely the frequencies, $\omega$, at which the atoms vibrate. The relationship between $\omega$ and $\{\nu,\mathbf{q}\}$ is called the __phonon band structure__, __dispersion relation__, or __dispersion curve__. It is commonly denoted as $\omega_\nu(\mathbf{q})$. It is very important because it summarises the relationship between the spatial and temporal properties of all phonons in one plot.  

More important quantities can be obtained from the phonon band structure. For example, the energy of the phonon, is simply given by  
$$
E_\nu(\mathbf{q})=\hbar\omega_{\nu}(\mathbf{q}).
$$
More examples related to the thermodynamical properties will be discussed in the next lab. 

### Dynamical matrix as a function of $\mathbf{q}$
Calculating the phonon band structure requires the dynamical matrix too, just like how we needed it for the finding the normal modes of molecules. In particular, we need to calculate the dynamical matrix for every wavevector $\mathbf{q}$. We will not go into the details. All you need to know here is that the dynamical matrix is going to be different for phonons with different wavevectors, and Quantum Espresso will calculate them for us.  


## Vibrations in crystals: Quantum Espresso (overview)
We will now learn how to calculate phonon band structures using Quantum Espresso. The overall procedure is similar to calculating the vibrational modes of molecules, with a few additional steps. 

The additional steps are required for reducing the computational cost. This calculation is very expensive in principle. This is because we need to calculate a new dynamical matrix $\mathbf{D}(\mathbf{q})$ for every phonon wavevector $\mathbf{q}$ for the entire phonon band structure. Quantum Espresso works around this problem smartly. We will not go into how it does it here. Essentially, it allows us to calculate the dynamical matrices for a smaller set of wavevectors, and use that to calculate the phonon band structure over a denser grid of wavevectors. 

Before going through the actual procedure, there are two things we need to tell you. Firstly, the workaround requires the use of two additional modules, `q2r.x` and `matdyn.x`. Secondly, there are two $\mathbf{q}$-grids that we need to specify. The first $\mathbf{q}$-grid is the one for which we calculate the dynamical matrices for. This grid is coarse. The second $\mathbf{q}$-grid is the one over which we calculate the phonon band structure. This grid is dense. Let's now begin our overview of the procedure.

The calculation has five stages. 

1. Perform a self-consistent calculation of the density and wavefunction. The module for this is `pw.x` as usual. 
2. Calculate the dynamical matrix on a coarser grid of wavevectors. The module for this is `ph.x`, which is the one we have used for molecules. 
3. Obtain a set of real space force constant matrices $\mathbf{K}(\mathbf{R})$ based on the dynamical matrices $\mathbf{D}(\mathbf{q})$ obtained in the last step. The module for this is `q2r.x`.  
4. Obtain the dynamical matrix over a denser grid of wavevectors. The module for this is `matdyn.x`.
5. Generate the phonon band structure plot. This will be done using Python. 

## Vibrations in crystals: Quantum Espresso (calculation)
We will now go through the calculations using carbon diamond as an example. 

### Step 1: run the `pw.x` calculation
!!! example "Task 4a - run `pw.x`"
    - Read the input file `01_CD_scf.in`. 
        - Again the variable `ecutrho` is set tighter. 
        - The prefix is defined explicitly as `'CD'`. This is useful when you run multiple calculations in different directories.
        - Other variables should be familiar to you.  
    - Run `pw.x` using this input file and save the output. Check that the job is completed. 
    - The other output files generated will be in `CD.save` since we set the prefix.

### Step 2: run the `ph.x` calculation 
!!! example "Task 4b - run `ph.x`"
    - Read the input file `02_CD_ph.in`.
    ```bash 
    phonons of Carbon diamond on a grid 
     &INPUTPH
      prefix = 'CD',         #(1)!
      asr = .true.           #(2)! 
      ldisp = .true.         #(3)!
      nq1=4, nq2=4, nq3=4    #(4)!
      /
    ```
        1. We have specified the same prefix as in the scf input file.
        2. Switch on the acoustic sum rule.
        3. `ldisp = .true.` tells Quantum Espresso to perform dynamical matrix calculation over a grid of q-points.
        4. `nq1`, `nq2` and `nq3` define our q-point grid. This is the coarser grid. 
    - This will produce a number of dynamical matrices. They have the same format as the `matdyn` file for the methane molecule. It will specify the $\mathbf{q}$-points at which the dynamical matrix is being evaluated within the file too. 
    - Note that one file could contain multiple $\mathbf{q}$-points because the dynamical matrices is the same at these $\mathbf{q}$-points due to symmetry. 

### Step 3: run the `q2r.x` calculation
Next we'll be generating the real space force constants from our grid of dynamical matrices using the `q2r.x` module.

!!! example "Task 4c - run `q2r.x`"
    !!! warning "Help file of `q2r.x`"
        `q2r.x` doesn't have a help file like the other codes. You need to inspect the source file in `PHonon/PH/q2r.f90` if you download a copy of the quantum espresso source code. The inputs are all described in a comment at the top of this file. On the mt-student server you can see this file at `/opt/build/quantum-espresso/q-e-qe-6.3/PHonon/PH/q2r.f90` [**replace this**].
    - Read the `q2r.x` input file. 
    ```bash
     &input
       fildyn = 'matdyn'     #(1)!
       zasr = 'simple'       #(2)!
       flfrc = 'CD444.fc'    #(3)! 
     /
    ```
        1. `fildyn` is the name of the FILe for writing the DYNamical matrix from `ph.x`. 
        2. `zasr` decides how the acoustic sum rule is enforced. If the ASR is not enforced, the lowest frequency phonons at the $\Gamma$-point will not have zero frequencies. This is not a physical behaviour, as we will explain later.    
        3. `filefrc` is the name of the FILE for writing the FoRCe constants in real space.
    
    - Run `q2r.x` with this input file now and save the output. This will run
      almost instantly.
    - This will produce the output file `CD444.fc` has the force constants for each pair of atoms in a 4x4x4 supercell. You do not need to know how to read this file now.


### Step 4: run the `matdyn.x` calculation
Now we want to use this to generate our normal mode dispersion. We'll be doing this with the `matdyn.x` code.

!!! example "Task 4d - run `matdyn.x`"
    !!! warning "Help file of `matdyn.x`"
        As with `q2r.x` this doesn't have a doc file describing its input variables. But you can check the comments at the top of its source file to get their details. On the mt-student server this is at `/opt/build/quantum-espresso/q-e-qe-6.3/PHonon/PH/matdyn.f90` [**change this**].
    - Read the input file `04_CD_matdyn-bands.in`.
    ```bash 
     &input
        asr = 'simple'            #(1)!
        flfrc = 'CD444.fc'        #(2)!
        flfrq = 'CD-bands.freq'   #(3)!
        dos=.false.               #(4)!
        q_in_band_form=.true.     #(5)!
     /
    8   #(6)!
     0.000 0.000 0.000 30   #(7)! 
     0.375 0.375 0.750 10
     0.500 0.500 1.000 30
     1.000 1.000 1.000 30
     0.500 0.500 0.500 30
     0.000 0.500 0.500 30
     0.250 0.500 0.750 30
     0.500 0.500 0.500 0
    ```
        1. `asr` to again tell it to try to force the acoustic modes at gamma to go to zero.
        2. `flfrc` to  give it the name of the file with the real space force constants from the `q2r.x` calculation.
        3. `flfrq` to give it the name of the output file to store its calculated frequencies.
        4. `dos=.false.` to tell it we're not calculating a density of states 
        5. `q_in_band_form=.true.` to tell it we want to calculate bands between high-symmetry points.
        6. The number of high-symmetry points which marks the high-symmetry path for calculating the phonon band structure. 
        7. The list of high-symmetry poinst with the number of points to calculate along each line, in the same way as we did for the electronic band structure.
    
    
    - Run `matdyn.x` using this input file now and save the output. Again this
      is very fast.
    - There's very little actual output from the code itself, but it will generate
      the files `CD-bands.freq` and `CD-bands.freq.gp`. Both of which contain
      the frequencies along the lines we requested.
    
### Step 5: plotting the phonon band structure 
Finally, we want to generate a plot of these frequencies. We could do that directly with the previous output: `CD-bands.freq.gp` is a multicolumn space separated list of the frequencies in cm-1.

!!! example "Task 5 - plotting the phonon band structure"
    Once you've done this, a Python script `VISUAL_PHONON_BAND.py` has been provided for generating the phonon band structure. If you are interested in how the plotting is done, you can read the script. 

    ??? success "What does the phonon band structure look like?"    

    ??? success "How many normal modes are there at each q-point?"


    <!--??? success "Inspect the phonon band structure along $\Gamma-L$, which two phonon bands become degenerate?"-->
    <!--    The two lowest energy phonon bands. Degenerate bands tend to appear in crystals with some symmetry. -->

    <!--??? success "Near which q-point do the three lowest energy phonon bands become approximately linear?" 
        The $\Gamma$ point. In 3D crystals, the three lowest energy phonon bands always become linear near the $\Gamma$ point. If this is not the case, you might need to increase the PW cutoff or k(q)-point sampling in `pw.x` and `ph.x`. -->

    <!--??? tips "Extra notes: use of plotband.x (for the interested)"
        It can be easier to use the `plotband.x` tool to help generate a plot. The steps are as follows:  
        1. Call this with `plotband.x < CD-bands.freq`.  
        2. This will then ask you for an Emin and Emax value - you should pick values equal to or below and above the numbers it suggests.  
        3. Then it will ask you for an output file name. Pick "CD-bands-gpl.dat" here.  
        4. You can then cancel further running of this code with `ctrl+c`. Note it has output the location of the high-symmetry points along its x-axis. -->

!!! example "Task 6 - visualising the phonons"
    - Go to the website https://interactivephonon.materialscloud.io/. This is a website for visualising phonons.


??? tip "Extra notes - what trick is being played by Quantum Espresso? (for those interested in the math)"
    Quantum Espresso calculate the function $D(\mathbf{q})$ at a coarse set of $\mathbf{q}$-points first. using the inverse Fourier transform, the force constant matrix in real space $K(\mathbf{R})$ is obtained. This is done with `q2r.x`. Then $K(\mathbf{R})$ is Fourier transformed again into the $\mathbf{q}$-space so that we can evaluate what $D(\mathbf{q})$ is over a denser grid of $\mathbf{q}$-points. This is done with `matdyn.x`. 

    Essentially we evaluate the dynamical matrix accurately at a few q-points, then estimate the rest using Fourier transformations. So instead of needing to carry out the DFPT calculation over a dense grid of $\mathbf{q}$-points, which will be computationally expensive, we can get away with a fewer number of DFPT calculations.  


??? tip "Extra notes - what is the acoustic sum rule? (for those interested in the physics)"
    The acoustic sum rule (ACS) states that the frequencies of acoustic phonon modes must be zero. This is due to the translational symmetry of crystals. Since the acoustic phonon modes at the $\Gamma$ point simply corresponds to rigid translation of the entire crystal, there must not be any relative displacements between the atoms from one unit cell to another. Hence, the vibration frequencies have to be zero. This rule is implemented in many different ways. You can look it up if interested. 

    For molecular calculations, switching on the ACS guarantees the lower energy modes corresponds to translations and have zero oscillation frequency, as we have seen earlier in the lab. 


------------------------------------------------------------------------------

Summary
-------

In this lab we have seen

- For a molecule:
    - How to use the `pw.x` code to calculate a converged density and
      wavefunction and then use these with the `ph.x` code to calculate the
      vibrational modes using DFPT.
- For a crystal:
    - How to use the `pw.x` code to calculate a converged density and
      wavefunction.
    - How to use these with the `ph.x` code to calculate the phonon modes on
      a grid of wavevectors.
    - How to transform these to obtain real-space force constants using the
      `q2r.x` code.
    - How to use the real-space force constants to obtain a phonon mode
      bandstructure along lines between high-symmetry points in the
      Brillouin zone.
    - How to plot the bandstructure with gnuplot in a similar way to how we
      plot the electronic bandstructure in a previous lab.
