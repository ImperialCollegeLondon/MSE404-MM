# k-points convergence test
The script we are using is a slight modification of the script used in Lab 2. 
The only change we made is from varying the energy cut-off to varying the
k-point grid density.

0. Modify the `kpt_build.py` and put in the correct template input file using
   `C_diamond_base.in` as reference. If you have any trouble doing this, 
   checkout `kpt_build_model.py`.

1. Run `kpt_build.py` to generate input files with different k-point grid
   densities.

   ```bash 
   python kpt_build_model.py 
   ```

2. Run `run.sh` to run the calculations with different k-point grid densities
   (remember to load the quantum-espresso module befor doing so), and to extract
   the total energy from the output files.

   ```bash
   bash run.sh
   ``` 
   The results are stored in `data_kpt.txt`.

3. Run `convergence_processing.py` to plot the total energy as a function of the
   k-point grid density and to find the converged values.

   ```bash
   python convergence_processing.py
   ```

