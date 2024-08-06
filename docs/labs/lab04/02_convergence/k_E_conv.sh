#!/bin/bash

template="C_diamond_base_kE.in"
repstr_k="xxxx"
repstr_E="eeee"

for val_k in {02..10..2}
do
for val_E in {20..100..20}
do
  echo "Running for k = $val_k and E = $val_E"
  inp="C_diamond_${val_k}_${val_E}.in"
  # We add the g here to replace every entry on the line.
  sed "s/$repstr_k/$val_k/g" $template > $inp
  sed -i "s/$repstr_E/$val_E/g" $inp
  pw.x < $inp &> ${inp%.*}.out_conv_kE
done
done

awk '/number of k points/{nkpt=$5}/kinetic-energy cutoff/{ekin=$4}
     /^!.*total/{print nkpt, ekin, $5}' *out_conv_kE > etot_v_nkpt_ekin.dat
