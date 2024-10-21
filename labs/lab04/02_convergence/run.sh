#!/bin/bash

# Run scf calculations.
for i in {001..015};
do
pw.x < scf.kpt.$i.in > scf.kpt.$i.out
done

nat=$(grep 'nat' *001.in | awk '{print $3}' | tr -d ', =')
nat=$(echo "$nat.0" | bc)

# Loop through files
for i in {001..015}; do
    # Extract 'kpts' value from input files
    kpt_value=$(grep '0 0 0' scf.kpt.$i.in | awk '{print $1}')

    # Extract 'Total Energy' value from input files
    final_energy_value=$(grep '!' scf.kpt.$i.out | tail -1 | awk '{print $5}')

    # Multiply 'final_energy_value' by 13.6
    final_energy_value=$(echo "$final_energy_value * 13.6" | bc)

    final_energy_value_per_atom=$(echo "scale=10; $final_energy_value / $nat" | bc)

    # Append values to output file.
    ######################################################################
    #                          Final energy in ev!                       #
    ######################################################################
    echo -e "$kpt_value $final_energy_value_per_atom" >> data_kpt.txt
done
