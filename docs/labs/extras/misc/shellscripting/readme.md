Bash and Shell Scripting
========================

Bash Shell
----------

There are often several different shells installed on a Linux system. `bash` is
probably the most common, and is typically the default for new users on any
system. `bash` is running in the terminal, and it interprets and executes the
commands we type in.

### Environment

#### User Configuration

On startup, such as when a new terminal is opened, or you connect to a remote
system over the network, `bash` will read several configuration files depending
on how it detects it is running.

- For login shells (shells where you are prompted to login when you start it),
  `bash` will read configuration options in `~/.bash_profile` or `~/.profile`
  if these exist. Note `~` is interpreted by the shell as the users home
  directory.
- For non-login shells (such as when you open the terminal emulator) `bash`
  will read configuration options in `~/.bashrc`.
- A lot of the time, for simplicity, many people will put all their
  configuration options in `~/.bashrc` and leave their `~/.profile` empty
  except for a command to read the `~/.bashrc` file.

These files can be used to configure many aspects of the shell, such as how the
prompt looks in the terminal.

#### Environment Variables

Many important aspects of how `bash` behaves are governed by environment
variables. These can be modified to your preference in your configuration
files. To see the current value of one of these variables you can do, e.g.
`echo $PATH`. A `$` symbol is used when referring to the data stored in a
variable, and `echo` is a built-in bash command that outputs something to the
terminal.

- `PATH`: this tells bash which directories to look in for executables so they
  can be run without typing in their full path. For example, we can just type
  in `gedit` in the terminal to launch the program, as the `PATH` variable
  contains the directory holding an executable of this name. Different
  directories are separated by ":".
    - To add the directory `~/.bin` to your path, you can add the line
      `export PATH="~/bin:$PATH"` to the `.bashrc` file. Here the `export`
      statement is used to allow the variable to be usable by any child
      processes rather than local to the script. Note when we are naming the
      variable we want to assign a value to we don't use a `$`. Also we just
      want to add a directory to the existing `PATH` so we add `:$PATH` to keep
      all the existing directories.

### Scripting

Shell scripts allow you easily string together sets of commands in a useful
way. For example, you could write a script that runs a calculation, parses the
important results to another file, and then generates a plot.

This course will only be able to briefly cover scripting. If you're
interested in developing more advanced scripts, I strongly suggest referring
to the
[Advanced Bash-Scripting Guide](http://tldp.org/LDP/abs/html/index.html).

#### The basics

A simple example script to output "Hello World!" is as follows:

```bash
#!/bin/bash
# Simple example script that outputs "Hello World!"

echo "Hello World!"
```

- The first line tells the system what command is used to interpret the
  contents of the script. This line should always begin with `#!` and then the
  path to the executable. For `bash` the executable will almost always be
  `/bin/bash`.
- Comments being with `#`.
- Commands can be entered otherwise as you would enter them to the terminal.
- Try creating a file called `hello.sh` with the contents listed above.
- You will need to make the script executable with `chmod u+x hello.sh`.
- Then you can run it by typing `./hello.sh`.

#### Variables

Variables don't need to be declared in any way, you can assign a value and
start using them. For example:

```bash
#!/bin/bash
# Example using variables.

var1="Hello"
var2="World!"

echo $var1 $var2
# Note the use of the $ symbol when we want to use the value stored in the
# variable.

# Any command can use these variables
dirname="tmpdir"
mkdir $dirname
```

You can also read a value from stdin using the `read` command as follows:

```bash
#!/bin/bash
# Example showing how the read command is used.

echo "Please enter some text:"
read user_text

echo "The text you entered was:" $user_text
```

##### Capitalization

Note: variable names are case sensitive in bash. The usual convention used is
that environment variables (such as `PATH`), and internal shell variables
(such as `BASH_VERSION`) are capitalized, while other variables are in
lowercase. Adopting this convention will ensure that you don't accidentally
overwrite any important variables in your scripts.

#### Command Substitution

Command substitution allows the result of a command to replace the command
itself, acting much like a variable in practice. For example:

```bash
#!/bin/bash
# Example of command substitution.

# This can be done by enclosing the command in $( )
echo "It is now $(date)."
# Note: the "date" command outputs the current date and time.

# This can also be done by enclosing the command in backtics ` `
echo "The files in this directory are:" `ls`
```

#### Conditional Statements

`if` statements can be used in bash, and many types of tests are possible.
You can test if the value stored in a variable equals something as follows:

```bash
#!/bin/bash
# Example using if statatements

echo "Please enter a yes or no: "
read user_response

# Note the spaces following `[` and before `]` are important.
if [ $user_response = yes ]
then
    echo "You entered yes."
elif [ $user_response = no ]
# "elif" is the same as "else if"
then
    echo "You entered no."
else
    echo "You didn't enter yes or no."
fi
# "if" statements are always ended with "fi".
```

We can also check, e.g., if a file or directory exists:

```bash
#!/bin/bash
# Check if the directory "tmpdir" exists, and if not, create it, then check
# if the file "tmpdir/testfile" exists, and if not, create it.

dirname="tmpdir"
filename="testfile"

if [ ! -d "$dirname" ]
# The "!" is logical negation
# -d tests that a file exists and is a directory.
then
    mkdir $dirname
fi

if [ ! -f "$dirname/$filename" ]
# -d tests that a file exists and is a regular file (i.e. not a directory).
then
    touch "$dirname/$filename"
fi
```
Many other tests of this form are possible.

Note, bash can also use the `[[ ... ]]` test construction rather than
`[ ... ]`, with the former offering some more options and functioning more like
other programming languages, though it won't work in many other (non-bash, or
older bash version) shells.

#### Arithmetic and Testing

To test numeric values, the `(( ... ))` construction can be used. For example:

```bash
#!/bin/bash
# Example showing the use of the (( )) construction

var1=4
var2=5
var3=8

if (( var1 + var2 > var3 ))
# Note within the (( )) don't use $ symbols before variables.
then
    echo "$var1 + $var2 > $var3"
else
    echo "$var1 + $var2 <= $var3"
fi

# We can also use the construction to perform basic arithmetic
var4=$(( var1 * var3 ))
echo "$var1 * $var3 = $var4"
```

#### For Loops

For loops iterate a variable over a range of values. For example:
```bash
#!/bin/bash
# Simple for loop example.

# This construction loops over a space separated list of values. A variable
# whose contents contains spaces would work in the same way.
for i in 1 2 3
do
    echo "Iteration number $i"
done
```

A for loops can be used to apply a command repeatedly to a set of files. For
example:

```bash
#!/bin/bash
# Example showing a simple for loop over a list of arguments used to convert
# a set of data files with comma separted columns to tab separated columns.

for csvfile in "$(ls *.csv)"
do
    datfile="$(basename $csvfile .csv).dat"
    # basename is a tool to strip the suffix from a file. Here we use it
    # to construct a new filename with the .dat extenstion.
    sed 's/,/\t/g' $csvfile > $datfile
done
```

Bash also supports numeric loop ranges using the syntax
`{start..finish..increment)`. For example

```bash
#!/bin/bash
# Example of for loop with numeric increments

for i in {0..6..2}
do
    echo "The value of i is $i"
done
```

#### While Loops

Bash also supports while loops, where a set of commands are repeated
continuously while a given condition is true. For example:

```bash
#!/bin/bash
# Example of while loop in bash

counter=0

# Here -lt means less than
while [ $counter -lt 10 ]
do
  echo $counter
  counter=$((counter + 1))
done
```

It is often convenient to use `while` loops to operate on every line of a file.
For example:

```bash
#!/bin/bash
# Example of while loop reading lines from a file

linenumber=0

while read line
do
  linenumber=$((linenumber + 1))
  echo "$linenumber: $line"
  # This will prepend each line of a file with its line number.
done < test.txt
```

#### Arguments

We can also pass arguments to bash scripts from the command line. For example:

```bash
#!/bin/bash
# Example using command line arguments. Call this script with several
# arguments.

echo "The first argument is $1"
echo "The second argument is $2"
echo "The number of arguments is $#"
echo "The full list of arguments is $@"
```

We can loop over arguments using a for loop as follows:

```bash
#!/bin/bash
# Example using a for loop to iterate over command line arguments.

for arg in $@
do
  echo $arg
done
```

#### Functions

We can also define functions in bash as we would in other languages. These
allow sections of code to be reused as needed and can help make a script easier
to follow.

A simple example is as follows:
```bash
#!/bin/bash
# Simple example of a function in a bash script

# First we define a function called "hello", that will output "Hello World!"
# when called.
hello() {
  echo "Hello World!"
}

# To call the function, we just use its name.
# Note function names in the script will take precedence over executable names
# from the path in the script.
hello
```
