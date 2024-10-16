#!/bin/bash

template="CO_charged_m1.in"

# Strings to be replaced
boxstr="xxxx"
xy_str="yyyy"
z1_str="zzz1"
z2_str="zzz2"

# Equilibrium C-O distance in Angstroms
codist=1.128

rm -f ${template%.*}_etot_v_box.dat

for A in {10..20..2}
do
  inp="${template%.*}_${A}.in"
  # We save the output filename to a variable also.
  out="${inp%.*}.out"

  # Again we use bc to get the cartesian components for each input.
  # We place the C and O atoms such that they are at same distance 
  # from the centre of the supercell 
  xy=$(echo "$A/2" | bc -l)
  zC=$(echo "($A - $codist)/2" | bc -l)
  zO=$(echo "($A + $codist)/2" | bc -l)
  # Display info about current calculation
  echo "Starting SCF calculation with box=$A angstroms"

  # Create new input file using consecutive instances of sed to replace strings in template
  sed "s/$boxstr/$A/" $template | sed "s/$xy_str/$xy/g" | sed "s/$z1_str/$zC/" | sed "s/$z2_str/$zO/"   > $inp
  pw.x < $inp &> $out

  # Check if it is converged and print outcome
  if grep -q "^!.*total" $out; then
	  echo "Converged"
  else
	  echo "Not converged"
  fi
  # awk is inside the loop this time, and we are appending to the data file
  # after each calculation completes.
  awk -v box=$A '/^!.*total/{print box, $5}' $out >> ${template%.*}_etot_v_box.dat
done
