#!/bin/bash
 pw.x < 01_C_diamond_scf.in &> 01_C_diamond_scf.out
 pw.x < 02_C_diamond_nscf.in &> 02_C_diamond_nscf.out
 bands.x < 03_C_diamond_bands.in &> 03_C_diamond_bands.out
