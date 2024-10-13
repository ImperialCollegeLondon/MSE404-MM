#!/bin/bash

template="C_diamond_base.in"
repstr="xxxx"

for val in {02..10..1}
do
  inp="C_diamond_${val}.in"
  # We add the g here to replace every entry on the line.
  sed "s/$repstr/$val/g" $template > $inp
done

