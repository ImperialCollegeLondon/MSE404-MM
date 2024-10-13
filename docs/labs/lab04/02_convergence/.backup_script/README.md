# Using Bash script

0. copy files to working dir
```
cp bash/* .
```

1. generate input files:
```
bash ./k_conv.sh
```

2. run pw.x for each input file:
```
for val in {02..10..1}
do
pw.x < $inp &> ${inp%.*}.out_conv
done
```

3. extract total energy from output files:
```
bash get_data.sh
```

4. plot
```
gnuplot plot_k_conv.gplt
```

# Using python

0. copy files to working dir
```
cp python/* .
```

1. generate input files:
```
python k_conv.py
```

2. run pw.x for each input file:
```
for val in {02..10..1}
do
pw.x < $inp &> ${inp%.*}.out_conv
done
```

3. extract total energy from output files:
```
python get_data.py
```

4. plot
```
gnuplot plot_k_conv.py
```
