# Extra: High Symmetry Points

!!! note "High symmetry points"
    If you have a particular structure and you want to find out which are the
    important k-points, then [this
    website](https://www.materialscloud.org/work/tools/seekpath) is a useful tool.

Finding appropriate high symmetry points and their labels for a band structure
plot is beyond the scope of this course, but generally you need to find the
Brillouin zone for the system you're interested in, along with the names of
the high symmetry points, and how these should be represented in terms of your
reciprocal lattice vectors.

These can often be looked up in a table for a given structure, while keeping
in mind that even for the same structure you may come across papers where some
zone boundary high symmetry points will be labelled differently. You can
access these tables online at e.g. the [Library of Crystallographic
Prototypes](http://aflow.org/CrystalDatabase/) by typing the name of the
mineral into the search box. This will give you the space group number, which
you could then use with e.g. the [Bilbao Crystallographic
Server](http://www.cryst.ehu.es/) to find the reciprocal space coordinates of
the high symmetry points and their labels. Another good reference is [TU
Graz](http://lampx.tugraz.at/~hadley/ss1/bzones) which has nice interactive
visualizations of the most common types along with their labels.

For the diamond lattice example, we might do this as follows:

- The diamond lattice is FCC. We can find find the space group number from
  <http://aflow.org/CrystalDatabase/A_cF8_227_a.html> and see that it is
  number 227.
- We can enter this number at <http://www.cryst.ehu.es/cryst/get_kvec.html>
  and find several high symmetry points. Often, you'll want to pick the same
  path that was chosen in some previous work to ensure you can reproduce it
  correctly. In
  [`02_C_diamond_nscf.in`](../../../lab04/03_bandstructure/01_diamond/02_C_diamond_nscf.in)
  we have set up an input file for the path `Γ-K-X-Γ'-L-X-W-L` where `Γ'`
  indicates the gamma point in a different Brillouin zone. Note - these
  labels don't exactly match points given in the table. Many points have
  a number of equivalent positions on the Brillouin zone surface, and often
  different conventions can be used for different materials with the same
  structure.



