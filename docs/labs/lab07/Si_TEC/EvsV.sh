#!/bin/bash

template_scf="Si_scf.in"
template_ph="Si_ph.in"
template_q2r="Si_qr.in"
template_matdyn="Si_matdyn.in"

repstr_a="xxxx"
repstr_i="yyyy"

# The "for" construction we've used previously only handles integers.
# So we set an initial value for A and its delta value as variables.
a1=5.30
da=0.04

for i in {00..5..1}
do
  inp_scf="Si_scf_${i}.in"
  inp_ph="Si_ph_${i}.in"
  inp_q2r="Si_qr_${i}.in"
  inp_matdyn="Si_matdyn_${i}.in"
  # bc is a calculator that can calculate expressions in the shell. We
  # pass it our calculation with echo and save the result as a variable.
  a=$(echo "$a1 + $i * $da" | bc)
  sed -e "s/$repstr_a/$a/" -e "s/$repstr_i/$i/" $template_scf > $inp_scf
  sed -e "s/$repstr_a/$a/" -e "s/$repstr_i/$i/" $template_ph > $inp_ph
  sed -e "s/$repstr_a/$a/" -e "s/$repstr_i/$i/" $template_q2r > $inp_q2r
  sed -e "s/$repstr_a/$a/" -e "s/$repstr_i/$i/" $template_matdyn > $inp_matdyn
  pw.x < $inp_scf &> ${inp_scf%.*}.out
  ph.x < $inp_ph &> ${inp_ph%.*}.out
  q2r.x < $inp_q2r &> ${inp_q2r%.*}.out
  matdyn.x < $inp_matdyn &> ${inp_matdyn%.*}.out
done
