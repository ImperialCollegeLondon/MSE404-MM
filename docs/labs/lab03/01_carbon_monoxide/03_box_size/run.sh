#!/bin/bash

# Run pw.x for each input file sequentially
for i in {10..35..5}; 
do
	pw.x < CO_$i.in &> CO_$i.out 
done

