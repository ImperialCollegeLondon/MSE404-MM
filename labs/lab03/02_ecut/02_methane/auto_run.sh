#!/bin/bash

template="CH4_base.in"
repstr="xxxx"

for val in {10..50..5}
do
  inp="CH4_${val}.in"
  sed "s/$repstr/$val/" $template > $inp
  pw.x < $inp &> ${inp%.*}.out
done

awk '/kinetic-energy/{ecut=$4}
     /^!.*total/{print ecut, $5}' *out > etot_v_ecut.dat
