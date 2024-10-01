#!/bin/bash

template="C_diamond_base.in"
repstr="xxxx"

for val in {02..10..1}
do
  inp="C_diamond_${val}.in"
  # We add the g here to replace every entry on the line.
  sed "s/$repstr/$val/g" $template > $inp
  # run pw.x?
  # pw.x < $inp &> ${inp%.*}.out_conv
done

# Extract the total energy and number of k points from the output files?
# awk '/number of k points/{nkpt=$5}
#      /^!.*total/{print nkpt, $5}' *out_conv > etot_v_nkpt.dat
