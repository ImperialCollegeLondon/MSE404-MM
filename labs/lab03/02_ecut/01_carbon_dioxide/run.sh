#!/bin/bash

# Run pw.x for each input file sequentially
for i in {10..40..5}; 
do
        pw.x < CO2_$i.in &> CO2_$i.out 
done

