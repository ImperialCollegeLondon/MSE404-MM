# Using python

0. copy files to working dir
```
cp python/* .
```

1. run pw.x for each input file:
```
pw.x < 01_C_diamond_scf.in &> 01_C_diamond_scf.out
pw.x < 02_C_diamond_nscf.in &> 02_C_diamond_nscf.out
bands.x < 03_C_diamond_bands.in &> 03_C_diamond_bands.out
```

2. plot
```
python plotbands_shifted.py
```
