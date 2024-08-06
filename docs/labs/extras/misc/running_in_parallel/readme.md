Running in Parallel
===================

All modern DFT codes are capable of running in parallel, provided they have been
compiled to do so. While you may be familiar with a single application using
several threads when you start it so that it runs faster, the parallelisation
scheme used by many DFT codes allows them to run on a large number of different
machines simultaneously to complete a single DFT calculation.

In the quantum espresso package, this is achieved through the
use of [MPI](https://en.wikipedia.org/wiki/Message_Passing_Interface) in
which a number of copies of a program are started at the same time which
can then pass information among themselves, communicating over the network
or faster interfaces such as infiniband, or running several on the same
machine where they each use a single core (or possible more than one, as
each process could in principle use several threads also). So a calculation
that takes say 10 minutes without using any parallelization, would take
(slightly more than) around 5 minutes using two parallel processes.

Throughout the course material, we make no mention of running in parallel.
This is just for simplicity, since it's one less thing for you (and me) to
worry about in the labs. The calculations we have you do in the labs and
homework assignments are small enough that this isn't necessary, but if you're
interested in doing more serious DFT calculations such as for a MSc project
then you should certainly start running your calculations in parallel.

Getting a version of espresso that can run in parallel
------------------------------------------------------

If you're on one of the mt-student servers, there is a separate module for
quantum-espresso compiled with parallel features enabled called
`espresso-mpi`. To load this, you'll also need to have the `openmpi` module
loaded first. It is set to conflict with the `espresso` module, so it will
generate an error if you try to load it while you have that loaded; if so
you can first unload that with `module unload espresso`.

> To load the parallel module on an mt-student server type `module load gcc
  mkl openmpi espresso-mpi`

If you have installed a VM on your laptop for this course, and have installed
quantum-espresso from Ubuntu repositories (e.g. via `apt`), this version
already has parallel features enabled.

Running your calculation in parallel
------------------------------------

There are several ways to start a parallel calculation, and if you're
using some HPC service, you such check their documentation for their
recommended approach.

To start a parallel calculation on mt-student or your own VM you can do the
following:

- Use the `mpirun` command, which is used to start a program that has been
  compiled with MPI enabled communication (if you run it on a normal program
  it will simply start several copies of that program at the same time).
- To tell it how many processes you want to start, you give it the `-np` flag
  followed by a number, such as `-np 2` to start two parallel processes.
- Then give it your program and input and output as usual.

Say for example, we have an input for a silicon calculation for `pw.x`
called `Si.in` and we want to save the output in `Si.out`:

- For the serial (non-parallel) calculation we would write
  `pw.x < Si.in &> Si.out`
- For a parallel calculation, if we wanted to use two parallel processes,
  we would write `mpirun -np 2 pw.x < Si.in &> Si.out`.

The majority of the codes that come with the quantum espresso package can
run in parallel in this manner.

You should be aware however, that planewave DFT calculations don't scale
linearly. Your calculation will get faster to a certain point, after which
if you add more parallel processes you'll slow your calculation down. This
can vary depending on the system and type of calculation you're doing, but
usually you'll see a reduction up to around 50 processes depending on the
parallelisation scheme (see below) and system involved.

Types of parallelisation
------------------------

If we do the above and write `mpirun -np 2 pw.x < Si.in &> Si.out`, we accept
the default strategy for parallelising the calculation. Different DFT codes
use different defaults, which all have their own advantages and disadvantages.
Using the default for quantum espresso is generally pretty good. The
differences between the different schemes are discussed in detail in the
quantum espresso documentation for `pw.x` but in general a planewave DFT code
can be parallelised in the following ways (which can be used in combination
provided this has been implemented):

- Over sets of calculations - if the calculation you have asked for involves
  running several similar calculations automatically, you can break up these
  calculations between your parallel processes. Quantum espresso offers this
  functionality for some types of calculations (such as for phonons) and
  refers to these sets as *images*.
    - You can set how many of these are used for a parallel calculation with
      `-nimage` or `-ni`. If you run with 20 processors and specify `-ni 2`,
      each image will use 10 processors. The default is 1.
- Over k-points - this scheme needs very little communication between
  processes and so offers very good scaling, as each k-point can be treated as
  effectively a separate calculation where results are added together at the
  end. It doesn't do as much as other schemes to reduce the memory
  requirements of the calculation, and if your calculation doesn't use many
  k-points, or you're calculating a molecule this may be a limited approach.
    - You can set how many parallel groups of k-points your calculation uses
      with the `-npools` or `-nk` flags in quantum espresso. The default is
      1.
- Over bands - This can cut down the amount of memory used by each process but
  requires a bit more communication between processes.
    - You can set how groups of bands are used for parallelising your
      calculation with the `-nband` or `-nb` flags. The default is 1.
- Over planewaves (FFT planes) - the plane wave basis set can be distributed
  across parallel processes. Quantum espresso does this very efficiently and
  this is its default parallelisation scheme. It will distribute planes of 3D
  FFT points to the parallel processes.
    - This is always turned on. Whatever the number of parallel processes that
      are left over after you specify other options will be used in this
      manner.
- Over task groups - if you have more parallel processes than FFT grid planes,
  you can redistribute the FFTs into task groups that allows for more
  efficient processing.
    - You can set the number of task groups with the `-ntg` or `-nt` flags.
      The default is 1.

There is also an overview of the options for the various quantum espresso
packages at <https://www.quantum-espresso.org/Doc/user_guide/node18.html>
if you'd like more detail.
