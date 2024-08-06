Additional Material: The Hydrogen Atom and Electron Spin
========================================================

As an example of a system where spin is important, let's look again at
hydrogen. In the molecule we have two hydrogen atoms bonded together,
and two electrons in total, and we can do this with a standard calculation.

But if we instead want to find the energy of an isolated H atom accurately
it's a little more tricky. A hydrogen atom has one electron, which means if we
treat it with doubly-degenerate bands, we would have a single half occupied
band. To calculate it like this we can treat it as a metal, and use some small
smearing, so that we allow partial occupation of bands. This is equivalent to
assuming that we have half an electron in each spin state. If we however
restrict it to being in one spin state or another, we may find a slightly
different energy (coming primarily from differences in how the DFT exchange
term is calculated in each case, with the latter being more physical).

The directory [`01_H1_metal`](01_H1_metal) has an input file for a single
H atom using a small smearing, while the directory [`01_H1_spin`](01_H1_spin)
has the same calculation, but with no smearing, and we have used the
input variables `nspin` and `tot_magnetization` to enable a spin polarized calculation.

- Run the input files in these two folders.
- Compare the total energy obtained in each case. In which case is the
  energy lower?
- Compare the energies of the lowest energy calculated bands.
    - Enabling smearing for the "metal" calculation will automatically add
      extra bands, but only the lowest energy band will be occupied in this
      calculation.
    - In the spin polarized calculation if you check the output you will see
      two sections for the band energies. One listing the energies of the
      "spin up" bands, and the other listing "spin down" bands. Since we
      have said we want 1 spin up electron in the calculation, the spin up
      band will be occupied and should be lower in energy than the
      unoccupied spin down band.
- Slightly above the band energies in the output of the spin polarized
  calculation, you'll see that quantum espresso also outputs the magnetic
  moment associated with each atom in the system. And slightly below the final
  total energy output, it will list the total magnetization of the system in
  Bohr magnetons per cell. The  measured value for hydrogen is 1. How close
  are you here?

