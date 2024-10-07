#!/bin/bash

# Run scf calculations.
for i in {001..010};
do
mpiexec pw.x < scf.mol.$i.in > scf.mol.$i.out
done

# Loop through files
for i in {001..010}; do
    # Extract 'ecutwfc' value from input molecule files
    ecutwfc_value=$(grep 'ecutwfc' scf.mol.$i.in | awk '{print $3}' | tr -d ', =')

    # Extract 'Total Energy' value from input molecule files
    final_energy_value=$(grep '!' scf.mol.$i.out | awk '{print $5}')

    # Multiply 'final_energy_value' by 13.6
    final_energy_value=$(echo "$final_energy_value * 13.6" | bc)

    echo "$ecutwfc_value $final_energy_value" >> data.txt
done
