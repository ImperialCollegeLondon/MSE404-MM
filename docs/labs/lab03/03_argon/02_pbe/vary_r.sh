#!/bin/bash

template="Ar2_base.in"
repstr="xxxx"

# delete the file if it exists already
rm -f etot_v_r.dat

for val in 3.3 3.4 3.5 3.6 3.7 3.8 3.9 4.0 4.1 4.2 4.5 5.0 
do
  echo "Calculating Ar dimer with r=${val}"
  inp="Ar2_r${val}.in"
  sed "s/$repstr/$val/" $template > $inp
  pw.x < $inp &> ${inp%.*}.out
  echo -en "${val}\t" >> etot_v_r.dat
  awk '/^!.*total/{print ecut, $5}' ${inp%.*}.out >> etot_v_r.dat
done

