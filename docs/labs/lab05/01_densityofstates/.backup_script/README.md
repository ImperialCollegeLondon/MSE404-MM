# Using Bash script

0. copy files to working dir
```
cp bash/* .
```

1. run pw.x for each input file:
```
pw.x < 01_C_diamond_scf.in &> 01_C_diamond_scf.out
pw.x < 02_C_diamond_nscf.in &> 02_C_diamond_nscf.out
dos.x < 03_C_diamond_dos.in &> 03_C_diamond_dos.out
```

2. plot
```
gnuplot plotdos_shifted.gplt
```

# Using python

0. copy files to working dir
```
cp python/* .
```

1. run pw.x for each input file:
```
pw.x < 01_C_diamond_scf.in &> 01_C_diamond_scf.out
pw.x < 02_C_diamond_nscf.in &> 02_C_diamond_nscf.out
dos.x < 03_C_diamond_dos.in &> 03_C_diamond_dos.out
```

2. plot
```
python plotdoss_shifted.py
```
