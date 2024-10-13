#!/bin/bash

# Extract the total energy and number of k points from the output files?
awk '/number of k points/{nkpt=$5}
     /^!.*total/{print nkpt, $5}' *out_conv > etot_v_nkpt.dat
