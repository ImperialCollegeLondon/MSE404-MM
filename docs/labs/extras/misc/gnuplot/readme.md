Plotting with Gnuplot
=====================

Gnuplot is a command-line driven open-source plotting utility, with many
features such as fitting, and 3D plotting available. You can install it on
ubuntu systems by typing `sudo apt install gnuplot`. The homepage is
[gnuplot.sourceforge.net](http://gnuplot.sourceforge.net), and a detailed
manual for the latest release is [also
available](http://gnuplot.sourceforge.net/docs_5.0/gnuplot.pdf). Gnuplot is
also readily scriptable. This allows you, for example, to incorporate it into a
bash script to automatically produce a file containing a plot of your results
after your calculation has finished.

To open gnuplot, simply type `gnuplot` in a terminal. You will see some
information regarding the version of gnuplot that has started, and finally a
gnuplot prompt: `gnuplot> `. You can enter various commands here to generate
and save plots.

For example:
- `plot sin(x)`
    - This will plot the sin function. The x values will range from -10 to +10
      by default and the y range will be automatically chosen to be -1 to 1.
- `plot cos(x), x + 0.1*x**2`
    - This will plot the cosine function in addition to the function
      y=x-0.1*x^2 .
- `plot "plot/example.dat"`
    - Plots the values listed in the example file `plot/example.dat`.
- `set title "My Results"`
    - This sets a title for the plot.
- `set xrange [-1:1]`
    - Sets the range of the x-axis in the plot. `yrange` can be set similarly.
- `set xlabel "Position (pm)"`
    - Sets the label for the x-axis in the plot. `ylabel` can be set similarly.
- `replot`
    - After changing the plot by e.g. adding a title, it is necessary to redraw
      the output plot. The `replot` command repeats the last plot command.

#### Outputing to a file

To output a plot to for example a pdf file, you need to set the gnuplot
"terminal" appropriately (the terminal setting determines the type of output
generated by gnuplot), set an output filename, and redraw the plot. Typically
many different terminals are available which allow ouput to e.g. postscript,
png, gif formats.

For example, to save a default plot of a sin function to a pdf:

1. `set terminal pdf`
2. `set output "sin_plot.pdf"`
3. `plot sin(x)`

#### Fitting

We can also define and fit functions within gnuplot. For example, to fit a
quadratic to the example data in `plot/example.dat` we can do the following
(here I assume gnuplot has been started from within the `plot` directory):

- `f(x)=a+b*x+c*x**2`
    - This defines the function in terms of a set of parameters.
- `fit f(x) "example.dat" via a,b,c`
    - This will do a least squares fit, and output the final parameter values
      along with standard errors.
- Note if no initial values for the parameters are specified, gnuplot will
  start each at 1. You can specify initial values before running the `fit`
  command as e.g. `a=-1;b=-1;c=0.5`. It is particularly important to give good
  initial guesses when fitting non-polynomial functions.
- `plot "example.dat", f(x)`
    - This will generate a plot of the data points together with the fit curve.
      You can also use this to try to find good initial guesses for parameters
      manually when fitting more complex functions.
- A summary of the fit results is automatically saved in the file `fit.log`.

#### Scripting

One can create scripts as a list of gnuplot commands entered in the same way as
would be done manually. Then `gnuplot scriptname` will execute the script and
exit. An example script to perform a quadratic to the data in `example.dat` and
generate a pdf plot of the data compared with the fit is given in
`plot/example.gpl`:

```
f(x)=a+b*x+c*x**2
a=-1;b=-1;c=0.5;
fit f(x) "example.dat" via a,b,c
set title "Example Gnuplot Plot"
set xlabel "Position (Bohr)"
set ylabel "Energy (Hartree)"
set term pdf
set output "example-gp.pdf"
plot "example.dat" with lines title "Results", f(x) title "Quadratic fit"
```

Try entering the `plot` directory and running this as `gnuplot example.gpl`.
You will see information on the fit output directly to the terminal, and the
files `fit.log` and `example-gp.pdf` will be generated. You can view the pdf
with the `evince` document viewer application that is installed by default on
ubuntu systems: `evince example.pdf`.
