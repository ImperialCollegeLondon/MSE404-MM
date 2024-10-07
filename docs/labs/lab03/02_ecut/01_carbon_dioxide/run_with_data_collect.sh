#!/bin/bash

# Run pw.x for each input file sequentially
for i in {10..40..5}; 
do
        pw.x < CO2_$i.in &> CO2_$i.out 
done

# Loop through files
for i in {10..40}; do
    # Extract 'ecutwfc' value from input molecule files
    ecutwfc_value=$(grep 'ecutwfc' CO2_$i.in | awk '{print $3}' | tr -d ', =') # Plane wave cutoff (Ry)

    # Extract 'Total Energy' value from input molecule files
    final_energy_value=$(grep '!' CO2_$i.out | awk '{print $5}') # Total energy (Ry)

    # Multiply 'final_energy_value' by 13.6
    final_energy_value=$(echo "$final_energy_value * 13.6" | bc) # Converting final energy to eV

    # Append values to output file.
    echo "$ecutwfc_value $final_energy_value" >> data.txt
done
