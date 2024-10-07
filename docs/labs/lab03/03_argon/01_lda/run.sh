#!/bin/bash

# Run scf calculations.
for i in {001..020};
do
mpiexec pw.x < scf.mol.$i.in > scf.mol.$i.out
done

# Loop through files
for i in {001..020}; do
    # Extract distance from input molecule files
    distance=$(grep 'Ar' scf.mol.$i.in | tail -1 | awk '{printf "%.3f", $4}')

    # Extract 'Total Energy' value from input molecule files
    final_energy_value=$(grep '!' scf.mol.$i.out | awk '{print $5}')

    # Multiply 'final_energy_value' by 13.6
    final_energy_value=$(echo "$final_energy_value * 13.6" | bc)

    echo "$distance $final_energy_value" >> data.txt
done
