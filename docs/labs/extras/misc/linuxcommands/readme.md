More Useful Linux Commands
==========================

While we cover the basic file and navigation commands in [Lab
1](../../../lab01/readme.md), There are many more commands that are available as
standard on a Linux system that can make your life much easier. In this file we
briefly go through some of the more commonly used ones and list a few useful
options or potential use cases for each. We'll revisit many of these during the
course as they'll help us automate some aspects of running and extracting data
from our calculations.

### Some additional useful file commands

As you are all on the same system and file space is not unlimited, it is good
to be aware of how much space you're using. Large calculations can often take
more disk space than you might have expected, so it's important to be able to
track down large files that are no longer needed if space becomes an issue.
The first few commands here can help you with this, beyond using `ls -lh` to
see how big files are.

### df

`df` is used to report file system disk space usage for all mounted
partitions. This is a useful way to check how much space is free on the disks
and partitions used to store files by the system. It also tells you which
disks are are real local disks: these are the ones with labels like
"/dev/sda2" under the filesystem heading, while disks mounted across the
network will be labelled with the server address in this field. The "tmpfs"
label indicates a temporary filesystem, usually this is stored in RAM while
the system is running, but can be used through a directory in the file
hierarchy.

You might notice a remote server directory is mounted to `~/homedir`. This
your windows home directory. You can copy files and folders there if you would
like to make the available on your Imperial windows account.

#### Useful Options

- `df -h` will list usage numbers in a human-readable format, i.e. instead of
  10240, it would say "10M".
- `df -T` will list the filesystem type e.g. ext3/ext4 for disks using a
  typical Linux filestructure, and nfs or cifs for disks mounted remotely over
  the network. It is useful to be able to identify filesystems that are
  mounted over the network since these are usually slower to use than local
  disks if your calculation produces many files. In these cases, if your home
  directory is on a remote filesystem you should run your calculation using a
  local disk and copy the results to your home directory once it is completed.

### du

`du` is used to estimate file space usage. Typing `du filename` will tell you
how big `filename` is (usually in kB by default). `du dirname` will tell you
how big all the files in directory `dirname` and its subdirectories are, with a
total at the end.

#### Useful Options

- `du -s` will only print the total size used for the argument. This is
  typically used to find how much space a given directory is using.
- `du -h` will print sizes in human-readable format.

### quota

`quota` is used to check disk usage and limits. Many shared Linux systems,
such as HPC systems, and the server you are using for this course, will impose
a limit on the disk space used by any individual users. The `quota` command
will tell you what those limits are and how much you have used.

#### Useful Options

- `quota -s` will show limits and usage in a human-readable format, i.e.
  instead of 3145180 it would say 3072M.

### head and tail

`head` outputs the first few lines of a file to the terminal, while
`tail` outputs the last few lines of a file to the terminal.

- `head filename` outputs the first 10 lines of `filename`.
- `head -5 filename` outputs the first 5 lines of `filename`.
- `tail filename` outputs the last 10 lines of `filename`.
- `tail -5 filename` outputs the last 5 lines of `filename`.
- `tail -f filename` starts tail in "follow" mode, where tail will repeatedly
  check for new data written to the file, and output them to the terminal.
  This is useful for following the output from a running calculation for
  example.

### chmod

#### Permissions

If you take a look at the output of `ls -l` you'll see the first column has a
mix of letters (usually d, r, w, x, and hyphens).

- The first character in this column indicates the file: `d` for directories,
  and `-` for regular files.
- The following nine characters indicate the permissions. These are in sets of
  three where each set indicates the permission for a set of users.
    - The first set of three characters are the permissions for the user that
      owns the file (listed in the third column).
    - The second set of three characters are the permissions of other members
      of the group that own the file (listed in the fourth column).
          - The default on many Linux systems is to create a group of the same
            name as the username, that contains only that user.
          - It's also possible for users to be added to several groups
            by the system administrator. This is useful on shared systems
            where a certain set of users want to share access to some set of
            files, but without giving access to everyone.
          - You can see what groups you are in by typing `groups` in the
            terminal.
    - The third set of three characters are the permissions for all other
      users.
- Within each set of three characters:
    - The first character is "r" if users in that set have permission to read
      the file and "-" otherwise.
    - The second character is "w" if users in that set have permission to write
      to the file and "-" otherwise.
    - The third character is "x" if users in that set have permission to
      exectute to the file, i.e. run it as a program.

#### chmod command

`chmod` can be used to change file permissions. This command can be invoked in
two different ways:

One can change the permission for a given set of users granularly:
- `chmod u+x filename` grants the user who owns the file execute permission.
  This is one of the main things you will be using `chmod` for as when you
  create a script in a text editor it will not be executable by default.
- `chmod g+rw filename` grants group members read and write permission.
- `chmod o-r filename` revokes other users read permission.
- `chmod a+x filename` grants all users execute permission.

One can use a set of three numbers to set the full list of permissions at once.
What each number corresponds to is listed in following table:

|  #  | Permission              | rwx |
|:---:|:------------------------|:---:|
|  7  | read, write and execute | rwx |
|  6  | read and write          | rw- |
|  5  | read and execute        | r-x |
|  4  | read only               | r-- |
|  3  | write and execute       | -wx |
|  2  | write only              | -w- |
|  1  | execute only            | --x |
|  0  | none                    | --- |


- To set the permissions of a directory so only the owner can access it in any
  way you could use `chmod 700 directoryname`.
- To set a script you have created so that others can use and execute it you
  could use `chmod 755 scriptname`.

### wc

`wc filename` will output the newline, word and byte counts for `filename`.

#### Useful Options

- `-l` to output the number of lines in the file.
- `-w` to output the word count for the file.

### grep

`grep` will print lines from a file or stdin which match a given pattern.

- `grep searchtext filename` will output all lines in `filename` which contain
  the text `searchtext`.
- `history | grep less` will output all lines in the command history containing
  `less`. This is useful for those times when you entered a complex command
  some time ago that you want to repeat it.

While `searchtext` in the first example above could be a particular word you
want to find an exact match of, `grep` will also interpret this as a *regular
expression* by default. This is somewhat similar to the wildcards you can use
in the terminal, but has a slightly different syntax and allows for much more
complex patterns.

#### Regular Expressions

This is a very deep topic, so we'll only cover a few of the more simple
examples. `man grep` has significantly more detail. The most useful symbols are
probably:

- `.` matches any single character.
- `*` the preceding item will be matched zero or more times.
- These are quite useful when combined to form `.*`, which acts in the same way
  as the terminal wildcard expression `*`.
- Note: to match the actual `.` or `*` symbols, you can escape them as `\.` and
  `\*`.

For example `grep "doc.*\.pdf" dirfiles.dat` will output all lines containing
strings that begin with `doc` and end with `.pdf`.

Note regular expressions can also be used in `less` (and hence `man`) when
searching for text with `/`.

#### Useful Options

- `grep -3 searchtext filename` will additionally output 3 lines before and
  after any lines containing `searchtext`. Any number of lines can be used
  here.
- `grep -v searchtext filename` will output all lines *except* those containing
  `searchtext`.
- `grep -r searchtext` will recursively search all files and folders starting
  from the current directory.

### cut

`cut` prints selected parts from each line a file. Mostly this is used with the
`-f` option which tells it to print only certain fields, and the `-d` option
which allows you to set the delimiter for fields (TAB is the default). It is
often useful to pipe (`|`) the output of `grep` into `cut` to parse data from
an output file.

For example, `cut -d ' ' -f 1 filename` will print the first word (separated by
spaces) on each line of `filename`, and `cut -d ',' -f '3-5' data.csv` would
print the 3rd, 4th and 5th columns of a csv data file.

#### Useful Options

- `-s` tells `cut` not to output lines which do not contain delimiters. This is
  useful for example if you have empty lines in the file you are parsing that
  you want to suppress.
- `--complement` will output the complement of the selected fields. For
  example, if we had a 5 column csv data file, `cut -d ',' -f '2-4'
  --complement` would output the 1st and 5th columns.

### awk

`awk` is a pattern scanning and processing language. It is typically used to
parse files in a similar manner to using `grep` combined with `cut`. `awk` is
very powerful but we will only cover some very basic operations.

- `awk '/regexp/{print}' filename` will output all lines in `filename`
  containing `regexp`. As with `grep`, regular expressions can be used in
  `regexp`.
- `awk '/regexp/{print $1" "$3}' filename` will output the first and third
  words in all lines containing `regexp`. Note by default `awk` uses spaces as
  the field delimiter.
- `awk 'BEGIN{i=0} /regexp/{i=i+1} END{print i}' filename` will output the
  number of lines in `filename` containing `regexp`.
- `awk '/searchtext/{printf "%f %f\n",$2-13.0,$4*10.0}' filename` will output
  for each line containing `searchtext`, the second field with 13.0 subtracted,
  and the fourth field times 10.

Hopefully these examples give you an idea of what is possible with `awk`. More
details and examples can be found with `man awk`.

#### Useful Options

- `-F` allows you to set the field separator. For example `awk -F','` would be
  useful for parsing a csv file.
- `-f program-file` tells `awk` to run the commands listed in `program-file`.
  This is useful if you have a complicated script so you don't need to type it
  all in directly to the terminal.

### sed

`sed` stands for stream editor. It allows you to perform basic text
transformations on a file (or from stdin). `sed` is very powerful and we will
only cover some very simple examples.

- `sed 's/regexp/replacement/' filename > newfile` will replace the first match
  of `regexp` on each line of `filename` with `replacement`. Here the first `s`
  stands for substitute, and is probably the most useful sed command. Note that
  `sed` outputs to stdout by default, so you should redirect this to a *new*
  file to save. **Do not try to redirect output to the same file you are
  reading**, as `>` will blank the file before `sed` can read anything from it.
    - `sed 's/^...//' filename > newfile` will remove the first three
      characters from every line of `filename`. Note `^` is used to match the
      beginning of a line.
- `sed 's/regexp/replacement/g' filename > newfile` will replace *every* match
  of `regexp` on each line of `filename` with `replacement`.
    - `sed 's/,/\t/g' data.csv > data.dat` would replace all commas with tabs
      (`\t` is a tab) in `data.csv` and save it in `data.dat`.
- The `-i` flag can be used to modify a file in-place.
  `sed -i 's/regexp/replacement/g' filename` will replace every match of
  `regexp` on each line of `filename` with `replacement`. `filename` itself
  will be modified. You can also specify a suffix here such that a backup
  will be created using that suffix: e.g.
  `sed -i.bak 's/regexp/replacement/g' filename` will do the replacement in-place
  but first backup the original file to `filename.bak`.

See `man sed` for more information.

### tr

`tr` is used to translate or delete characters. It always reads from stdin, and
outputs to stdout. This means that to use it with a file, we need to redirect
the file to stdin using `<`.

- `tr 1 2 < test.dat` would output the contents of `test.dat` with all the 1s
  replaced by 2s.
- `tr abc ABC < test.txt` would output the contents of `test.txt` with any 'a'
  replaced by 'A', 'b' by 'B' and 'c' by 'C'.

It also accepts some special input such as

- `[:space:]` to match whitespace (both single and continuous).
- `[:punct:]` to match punctuation
- `[:lower:]` to match lower case letters
- `[:upper:]` to match upper case letters

For example:

- `tr [:lower:] [:upper:] < test.txt > test_upper.txt` would to create a new
  version `test.txt` converted to uppercase.
- `tr [:space:] '\n' < test.txt` would convert all spaces to newlines.

#### Useful options

- `-d` deletes matching characters. For example, to output a file with all
  punctuation removed we could do `tr -d [:punct:] < test.txt`

### find

`find` is used to search for files in a directory hierarchy. Most commonly this
is used with the `-name` option to search for files with a particular name.
Wildcards can be used in the search. Note: the first argument to `find` should
be the path to search, e.g. `find /etc` to search for files in the `/etc`
directory or `find .` to search for files in the current directory.

- `find . -name "*.cpp"` will find all files ending in `.cpp` in the current
  directory (`.`) and its subdirectories.

See `man find` for more information.

### diff

`diff` is used to compare two files. This is useful if for example, you want to
see what changes have been made in a new version of a file.

- `diff file1 file2` will output the lines which differ between the two files.
  The lines from `file1` will be prepended with `< ` and the lines from `file2`
  with `> `.

#### Useful Options

- `-q` will report only whether the two files differ and will not output the
  differences.
- `-r` will recursively compare files in subdirectories.
- `-y` will output the two files side by side in two columns.
- `-W` will allow you to set how wide the output is (130 columns by default).
  This is particularly useful with the `-y` option.

### sort

`sort` is used to sort lines of text files. For example, if we had a file
called `users.txt` which contained a list of names, then `sort users.txt` would
output (to stdout) the list sorted alphabetically. This is often useful
combined with other commands. For example to generate a sorted list of all
words in a file you can do `sed 's/ /\n/g' filename | sort`. Here the `sed`
command replaces all spaces with new lines, so we have one word per line, and
then we use this as input to the sort command.

#### Useful Options

- `-n` will sort numerically rather than alphabetically. For example,
  `du -s * | sort -n` will generate a listing of files and directories sorted
  by size.
- `-h` will use a human numeric sort, allowing numbers such as 2K and 1G to be
  sorted. For example, `du -sh * | sort -h` will generate a listing of files
  and directories sorted by size, but in human readable format.
- `-u` will output only the first of an equal run. For example
  `sed 's/ /\n/g' filename | sort -u` will generate a sorted list of all words
  in `filename` with each listed only once.
- `-f` will till `sort` to fold lower case to upper case characters.
- `-r` will output in reverse order. This is useful for numeric sorts where you
  often want to have the largest numbers at the top.

### uniq

`uniq` is used to report or omit repeated lines in a file. By default it will
take any file or input from stdin, and output it with duplicated lines omitted.
For example if we had a text file `test.txt` with

```text
a
a
b
b
b
c
d
d
```

Running `uniq test.txt` would output

```text
a
b
c
d
```

#### Useful Options

- `-c` prefixes each line of output with a count of its number of occurrences.
  We could for example, take the `sort` example to sort all words in a file,
  and expand it to generate a word count of all words in a file:
  `sed 's/ /\n/g' filename | sort | uniq -c`. We could add `| sort -n` to the
  end of this to sort words in order of frequency, and we could use `tr` before
  this to remove all punctuation.
- `-i` tells `uniq` to ignore differences in case when comparing.

### tar

`tar` is an archiving utility used to create an archive of files, i.e. generate
a file containing many other files. This is usually used to create compressed
bundles of files on Linux, in a similar way to zip file archives (note zip and
unzip are usually available on Linux also, but compressed tar archives are more
commonly used).

#### Creating Archives

The `-c` flag indicates to `tar` that you want to create an new archive.

- `tar -cvf archive.tar file1 file2 dir1` will create an (uncompressed) archive
  called `archive.tar` of the named set of files or directories. Here:
    - `-v` is for verbose mode - it will list the files which are added. This
      is not really necessary, but is useful so you can be sure you are adding
      the files you intended to.
    - `-f` is used to specify the archive name. Here we have called it
      `archive.tar`.
- `tar -czvf archive.tar.gz file1 file2 dir1` uses the additional `-z` flag to
  compress the archive using `gzip` compression. The extension `.tar.gz` or
  `.tgz` is typically used to indicate this type of file. Several other
  possible compression algorithms could be used instead:
    - `-j` will use `bzip2` compression. This typically results in slightly
      smaller file size than `gzip`, but can take slightly longer. Files
      created like this usually use the extension `.tar.bz` or `.tbz` or some
      other similar variation.
    - `-J` will use `xz` compression. This is a very effective compression
      algorithm, resulting in very small files, particularly for archives
      containing a lot of text. This option can take quite a lot longer than
      gzip for large archives however. Files created with this option usually
      use the extension `.tar.xz` or `.txz`.

#### Extracting Archives

The `-x` flag indicates to `tar` you want to unpack an archive.

- `tar -xvf archivename` will uncompress and unpack a tar archive called
  `archivename`, automatically detecting what kind of compression was used (if
  any). Again the `-v` isn't necessary, but is useful.

#### Listing Archive Content

The `-t` flag will tell `tar` to list the contents of an archive.

- `tar -tvf archivename` will list the contents of the tar archive
  `archivename`, again automatically detecting what kind of compression was
  used.

### bc

`bc` is an arbitrary precision calculator language. This can be used from
the command line by piping expressions to it. For example:

- `echo "2+2" | bc`
- For floating point operations, you should set the scale which defines how
  many digits following the decimal points are output:
  - `echo "scale=10; 1.412*27.211" | bc`
  - `echo "scale=10; sqrt(2)" | bc`
- You can tell bc to load the standard math library with the `-l` flag. This
  will also set the scale to 20 by default. This makes several additional
  functions available such as `s(x)`, `c(x)` and `a(x)` for the sine, cosine
  and arctan in radians. So `echo "4*a(1)" | bc -l` will output the first 20
  decimal places of pi.


