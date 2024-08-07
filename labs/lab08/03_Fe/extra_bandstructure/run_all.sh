# run SCF calculation
pw.x < Fe.spin.scf.in > Fe.spin.scf.out
# run NSCF calculation
pw.x < Fe.spin.nscf.in > Fe.spin.nscf.out
# run bands calculation to extract band informtaion 
# Note that we have to do this for each spin channel
#
# Spin Channel Up
bands.x < Fe.spin.bands.up.in > Fe.spin.bands.up.out
# Save data to bands.out.up.gnu
cp bands.out.gnu bands.out.up.gnu
# Spin Channel Down
bands.x < Fe.spin.bands.dn.in > Fe.spin.bands.dn.out
# Save data to bands.out.dn.gnu
cp bands.out.gnu bands.out.dn.gnu
# Now plot the bands. Output is Iron_bands.png
gnuplot plotbands_shifted.gplt
