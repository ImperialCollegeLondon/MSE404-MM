#!/bin/bash

# Run scf calculations.
for i in {1..10};
do
    pw.x < Ar2_$i.in > Ar2_$i.out
done
