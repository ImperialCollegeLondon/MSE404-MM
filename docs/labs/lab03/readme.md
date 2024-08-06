Convergence and Importance Parameters
=====================================

In this lab we'll continue looking at molecules, and we'll also be going
through how to define various input parameters and how to check how well
converged your calculations are.

**Reminder** Don't forget to copy the `lab03` folder from `/opt/Courses/MSE404/lab03`
to your home directory.

Linux Recap
-----------

Before you start, if you can't remember how to do something from the command
line, this [cheat
sheet](https://www.slideshare.net/NoFernndezPozo/unix-command-sheet2014) may
come in useful. Or you can always refer back to [lab 1](../lab01/readme.md).

Although it's not essential, at some point you may wish to transfer files between
different computers. You can do this using the `scp` command. This
[website](https://linuxize.com/post/how-to-use-scp-command-to-securely-transfer-files/)
explains how to use it. There are also some other useful commands which we
haven't discussed [`here`](../extras/misc/linuxcommands/readme.md).


Pseudopotentials
----------------

As discussed in lectures, pseudopotentials are used to approximate the core
potential. For each atomic species, you need a file which describes
the approximation you want the DFT code to use for that species. As
you've seen we set the input files to look in the current directory for the
required pseudopotential file. This is done by setting `pseudo_dir = '.'` in the
`CONTROL` section of the input, where `.` represents the current directory.
We then use a link to this file from some central directory rather than making
multiple copies of the pseudopotential for different calculations. This week
we'll be doing some different molecules. As always we need to make sure each
of the atomic species have a pseudopotential listed in the input file.

An alternative way to do this would be to specify the central pseudopotential
directory directly in the input file. We've set this in the input file in
[`01_carbon_dioxide/CO2.in`](01_carbon_dioxide/CO2.in).

- Take a look at the directory contents and you'll see there's only an
  input file there.
- Compare this to one of the input files from [lab 2](../lab02/readme.md) and
  run the calculation to check it works. Note that you do not need to edit the
  file.
- Take a look at the pseudopotential file we've used for oxygen. The header
  has some useful information regarding how the pseudopotential was generated,
  such as what states are included, and what approximations are used for
  exchange and correlation.


Plane-wave energy cut-off
-------------------------

Regardless of the type of system you're looking at, you'll need to check how
well converged your result is (whatever it is your calculating) with respect
to the plane-wave energy cut-off. This governs how many plane-waves are
used in the calculation.

- In Quantum Espresso this is controlled with the parameter `ecutwfc`.
- Different systems will converge differently - for example you shouldn't
  expect diamond and silicon to be converged to the same accuracy with the same
  energy cut-off despite having the same structure and same number of
  valence electrons.
- Different pseudopotentials for the same atomic species will also converge
  differently. Often pseudopotential files will suggest an energy cut-off.
- Different calculated parameters will converge differently.
    - You should be particularly careful when calculating parameters that
      depend on volume, as the number of plane-waves for a given energy
      cut-off is directly proportional to the volume so this can introduce
      an additional variation. We'll see more about this later.

An example showing the total energy convergence with respect to energy cut-off
is in the [02_ecut/01_carbon_dioxide](02_ecut/01_carbon_dioxide) directory. We
have already set up a series of input files which are all identical except we
systematically increase the value of `ecutwfc`.

Examine one of these input files. You'll notice we've specified some
additional variables:

- In the `CONTROL` section we've added the setting `disk_io = 'none'`. This
  suppresses the generation of the wavefunction file and the folder with the
  charge density etc. If we only care about the total energy we don't need to
  generate these files, and the calculation will run a little faster if it's
  not trying to write unnecessary files as disk I/O can often be a bottleneck
  in calculations.
- In the `ELECTRONS` section we have added the `conv_thr` setting, though it
  is set to its default value. This variable controls when the self
  consistency cycle finishes. You should be aware of this variable, as there
  is little point in trying to converge to greater accuracy than we are
  converging self-consistently.

So now we want to run `pw.x` for each of these input files. **Don't forget
you'll need to load the espresso module and its dependencies before the
`pw.x` command will work.**

It's a little tedious to type out the `pw.x` command for each input file in
turn manually. Instead we can automate this with a small script.

Shell scripting
---------------

A shell script is a collection of Linux shell commands in a file, and when
that file is executed (in the same way as any Linux command is executed) then
those commands are executed. We could make a script that explicitly runs each
input file with `pw.x` in turn as follows:

```bash
#!/bin/bash

# Run pw.x for each input file sequentially
pw.x < CO2_10.in &> CO2_10.out
pw.x < CO2_15.in &> CO2_15.out
pw.x < CO2_20.in &> CO2_20.out
pw.x < CO2_25.in &> CO2_25.out
pw.x < CO2_30.in &> CO2_30.out
pw.x < CO2_35.in &> CO2_35.out
pw.x < CO2_40.in &> CO2_40.out
```

If we save this in a file called say `run_set.sh` (e.g. by copying the above 
to whichever text editor you prefer such as `gedit`), the simplest way to run the
script is using `bash run_set.sh`, where bash is the default shell on the system 
we're using. This means that the same commands you type in a terminal can be used
in a bash script. In other words this script does the same thing as if we typed
each line directly. It is also possible to execute the script directly using
`./run_set.sh`, which would rquire you to set the file to executable, using
[`chmod`](../extras/misc/linuxcommands/readme.md#chmod).

There are also a number of features available in `bash` to make scripts more
general than an explicit set of commands. For example, let's say we want to
instead find whatever input files are in the current directory, and run
`pw.x` with them, saving the output in an appropriate file, we could make our
script as follows:

```bash
#!/bin/bash

# Loop over files in the current directory with extension .in
for input_file in $(ls *.in)
   # For each file run pw.x and save output in file with same 
   # root name as input file but with .out extension
   do
      pw.x < "$input_file" &> "${input_file%.*}.out"
done
```

In this script, we save the output of an `ls` command listing all the input
files with `.in` extension as a variable. Then we loop over those input files, 
running `pw.x` with each, and saving the output in a file that has the same name as the input
file but with the extension `.out` instead of `.in`. The part
`${input_file%.*}` returns the value stored in the `$input_file` variable
but with the extension stripped away. This lets us make our script much more
general: if we add additional input files, they'll automatically be picked
up and run without us needing to modify our script.

### _Task_

 - Save this second script as `run_all.sh` and run it within the directory
   with the carbon dioxide input files using `bash run_all.sh`.
 - Check one of your output files to make sure Quantum Espresso ran
   as expected. If your file contains only an error message, you likely forgot
   to load the modules required to use Quantum Espresso on our server.

Once the calculations have run, we want to see how the total energy changes
with the input plane-wave energy cut-off. We could go through each output
file and find the resulting total energy and gather them in a file, but again
that would be tedious so we can write a script to do it for us (or extend
our earlier script to also do this).

There are two different commands we could use to extract the resulting total
energy from file. The first is often simpler to use and is called
[`grep`](../extras/misc/linuxcommands/readme.md#grep). For example we could use
the following to print the line containing the final total energy from each
output file using the fact that `pw.x` helpfully starts this line with a `!`:
`grep '^!.*total energy' *out`. In particular, we are searching for a line which
has the symbol `!` at the beginning `^` and that also contains the string `total
energy` in all files whose names end with `out`. Try running this now in the
directory containing your output files.

The other command we could use is
[`awk`](../extras/misc/linuxcommands/readme.md#awk), which is more powerful, but
also more complicate to use, and lets us pick out both the total energy value
and the energy cut-off that was used with a single command. We could use this in
a simple script as follows:

```bash
#!/bin/bash

# In all files ending with 'out', find line with 'kinetic-energy' 
# and save the fourth column/word into variable 'ecut'. 
# Then in the same file also find a line beginning with '!' and that contains
# the string 'total' and print the value stored in 'ecut' and the
# fifth column/word from this line.
awk '/kinetic-energy/{ecut=$4} /^!.*total/{print ecut, $5}' *out
```

Here we use the `awk` command to find the line with `kinetic-energy` and save
the fourth word (`$4`) to a variable called `ecut`. Then when it finds a line
that starts with a `!` and has the word `total`, it will output the value
stored in the `ecut` variable, followed by the fifth word (`$5`) on that line.

### _Task_

 - Save this as a script called `etot_v_ecut.sh`, run it in the directory with
   the output files and save the output as `etot_v_ecut.dat`. You can do this
   by typing `bash etot_v_ecut.sh > etot_v_ecut.dat`.
 - Look at the output (e.g. using `gedit`) and you can see to how many significant
   figures the total energy is converged for a given energy cut-off.

We can make things even easier if rather than manually generating all the input
files, we modify the script to take a base input file and modify it for each
calculation, run it, and finally parse the output to a data file.

This is in the `02_ecut/02_methane` directory.

The script is as follows:
```bash
#!/bin/bash

# Original filename
template="CH4_base.in"
# String to be replaced
repstr="xxxx"

# Loop from 10 to 50 in steps of 5. These are the values of the energy cut-off
for val in {10..50..5}
do
  # Define input file name. Assign at each file a specific name
  inp="CH4_${val}.in"
  # Substitute "xxxx" string in original file with the current energy cut-off value
  # and paste the result into new input file
  sed "s/$repstr/$val/" $template > $inp
  # Run pw.x on the current input file
  pw.x < $inp &> ${inp%.*}.out
done

# Extract for each output file the values of the energy cut-off
# and the final total energy
# Paste results in etot_v_ecut.dat file
awk '/kinetic-energy/{ecut=$4}
     /^!.*total/{print ecut, $5}' *out > etot_v_ecut.dat
```

Here we've combined the two scripts we created above, and also automated the
generation of input files using the
[`sed`](../extras/misc/linuxcommands/readme.md#sed) command. This can be used to
search for and replace some text in a file. We have set up a template input file
[`CH4_base.in`](02_ecut/02_methane/CH4_base.in) where we have used the
placeholder text `xxxx` as the text we'll search for and replace with energy
cut-off we want for our input file. The bash construction `for val in
{10..50..5}` will create a loop where the value stored in the variable `$val`
runs from 10 to 50 in steps of 5.


### _Task_

- Run this script and see what input and output files are generated.



Plotting with Gnuplot
---------------------

It's often useful to be able to generate a quick plot when testing the
relation between variables. `gnuplot` is a useful tool for generating plots,
particularly as it is also scriptable, so we could extend our earlier script
to automatically generate a plot from from the extracted data. There is a more
detailed overview in the [gnuplot section](../extras/misc/gnuplot/readme.md).

We can launch gnuplot by typing `gnuplot` in a terminal. Once it opens, we
can, for example, plot a data file by typing `plot "etot_v_ecut.dat"`
(provided we are in the directory containing that file). This will only plot
the points by default. If we want to join these up with lines we can type
`plot "etot_v_ecut.dat" with linespoints`, or `p "etot_v_ecut.dat" w lp`.
If you want to specify different columns of data, you can do for example
`p "etot_v_ecut.dat" u 2:1 w lp` to plot column 2 vs. 1 instead of 1 vs. 2
(which is the default).

### _Task_

- Generate plots of the difference between the total energy and its most
  converged value versus plane wave energy cut-off for both methane and
  carbon dioxide.
    - Instead of working out energy differences yourself, you can plot the
      difference in energies directly in gnuplot, for example
      by doing `eref=-100; p "etot_v_ecut.dat" u 1:($2-eref) w lp`, where 
      in this case we are calculating the energy relative to -100.
      You should change the value of `eref` from -100 to the energy of
      the calculation with the highest plane-wave cut-off.
    - It can also be useful to plot convergence on a logscale, which you can do
      in gnuplot by typing `set logscale y` before you type the above command.
- How does the behaviour compare between the two molecules?


Exchange & Correlation Functional
---------------------------------

The exchange and correlation is a key part of DFT. The functional we use
determines how we approximate all the many body interactions. By default,
within Quantum Espresso, the exchange and correlation functionals that are
used are taken from the header of the pseudopotential file as mentioned above.
It's possible to override this using the `input_dft` variable in the system
section, but it's best to use the same approximation as was used in the
pseudopotential generation.

The exchange correlation can have a big impact on a number of parameters.
In this example we consider the example of an argon dimer.  By varying the bond
length between two argon atoms we can plot binding energy curves. This is 
another example where shell scripting can be very useful, since we want to run
several similar calculations. We can modify the script we used for the cut-off
energy convergence to do this. We again define a template input file, found in
[`Ar2_base.in`](03_argon/01_lda/Ar2_base.in) and then use the following script,
which can be found in the `03_argon/01_lda` directory. In this case we want more
data points in some places than others, so we directly specify each distance we
want to calculate. There are also some other differences from the previous script -
in this case we take the value of _r_ that we set and directly print it using the
command `echo`, rather than extracting it from the Quantum Espresso output, but
otherwise the same principles apply.

The script for LDA is as follows:
```bash
#!/bin/bash

# Original filename
template="Ar2_base.in"
# String to be replaced
repstr="xxxx"

# delete the file if it exists already (-f)
rm -f etot_v_r.dat

# Loop over bond lengths
for val in 3.0 3.1 3.2 3.3 3.4 3.5 3.6 3.7 3.8 4.0 4.2 4.5 5.0
do
  # Print message to output
  echo "Calculating Ar dimer with r=${val}"
  # Define input file with specific name
  inp="Ar2_r${val}.in"
  # Substitue 'xxxx' string with current bondlength value
  # and paste result in new input file
  sed "s/$repstr/$val/" $template > $inp
  # Run pw.x on the new input file and save results in 
  # output file with same rootname as input file
  pw.x < $inp &> ${inp%.*}.out
  # Save current bondlength value in etot_v_r.dat file
  echo -en "${val}\t" >> etot_v_r.dat
  # Extract cut-off and total energy and append to etot_v_r.dat file
  awk '/^!.*total/{print ecut, $5}' ${inp%.*}.out >> etot_v_r.dat
done
```

### _Task_

- Run this script and see what input and output files are generated.
    - In order to save computer time the box size and cut-off have been
      chosen for speed rather than accuracy. If you have time, try converging
      these values properly. On the other hand, if the calculations take too
      long then reduce the number of data points by removing some distances.
- See how the total energy varies with the distance between atoms.

Now we want to try a different functional.  In the directory `03_argon/02_pbe`
an input file is present that is identical to the previous calculation, except
we specify a different pseudopotential here: `Ar.pbe-n-rrkjus_psl.1.0.0.UPF`.

### _Task_

- Look at this file, and compare the header section to the pseudopotential you
  used previously. You'll notice a line mentioning "PBE Exchange-Correlation
  functional" in the pseudopotential we're using here, where it said "PZ
  Exchange-Correlation functional previously".
    - Here PZ refers to a particular parametrisation of the Local Density
      Approximation (LDA) - the simplest approximation for the exchange and
      correlation, while PBE is more advanced approximation (though not
      necessarily better performing), which is classed as a Generalized
      Gradient Approximation (GGA). The letters PZ and PBE are the initials of
      the authors of the papers in which the particular approximations were
      published.
- Run the script for PBE.
- Plot distance vs. length for both PBE and LDA using the gnuplot script
  `plot_ar.gp` in the `03_argon` directory, by typing `gnuplot plot_ar.gp`. This
  will generate a file called `ar_dimer.eps`. You can view this file for example
  using the program `evince` if you type `evince ar_dimer.eps`. How do the two
  functionals compare?

------------------------------------------------------------------------------

Summary
-------

- In this lab we looked at defining pseudopotentials, checking the convergence of
  the total energy with respect to the plane-wave energy cut-off, and the effect
  of exchange and correlation functional.
    - Convergence of any parameter is done by systematically varying the
      corresponding calculation parameter and looking at how the result changes.
- We saw how we can use bash scripts to automate this process.
    - This means we don't need to manually create a number of almost identical
      input files, and manually go through each one to find the values we
      want.
    - We can use a bash `for` loop to perform a calculation for a number of
      input files.
    - We can use `grep` or `awk` to parse results or parameters from our
      output files.
    - We can use `sed` to replace values in a template input file.
- We can quickly generate a plot of a data file with `gnuplot`.

------------------------------------------------------------------------------

Extra 
------

### _Box-Size_

For systems which are not periodic in three dimensions, such as molecules which 
we are calculating within a box, we also need to test the convergence with respect
to the box size. In other words, we need to check that the spurious interaction
between periodic images is sufficiently small for the accuracy we desire. We would
also have to do something similar with for example 2D systems, where we would have
empty space in one direction.

### _Optional Task_

- Create a folder called `04_spacing/01_methane_20` and test the convergence
  of total energy versus box dimension for an energy cut-off of 20 Ry. And
  then create a folder called `04_spacing/01_methane_60` and test the
  convergence of total energy versus box dimension for an energy cut-off of 60
  Ry.
    - How do these compare?


### _Advanced Shell scripting_

If you're interested in reading more about scripting, you can take a look at the
[shellscripting section](../extras/misc/shellscripting/readme.md). We'll provide
plenty of examples and keep things relatively simple in the course, but you may
find some of the more advanced functionality useful in your own work.
