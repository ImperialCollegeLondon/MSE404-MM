# k-points convergence test
The script we are using is a slight modification of the script used in Lab 2. 
The only change we made is from varying the energy cut-off to varying the
k-point grid density.

1. Run `kpt_build_model.py` to generate input files with different k-point grid
   densities.

   ```bash 
   python kpt_build_model.py 
   ```
2. Run `run.sh` to run the calculations with different k-point grid densities,
   and to extract the total energy from the output files.

   ```bash
   bash run.sh
   ``` 

   The results are stored in `data_kpt.txt`.
3. Run `convergence_processing.py` to plot the total energy as a function of the
   k-point grid density.

   ```bash
   python convergence_processing.py
   ```

