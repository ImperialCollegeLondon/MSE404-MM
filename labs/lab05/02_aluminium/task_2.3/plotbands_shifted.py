import numpy as np
import matplotlib.pyplot as plt

# Reading the data from the file
data = np.loadtxt("bands.out.gnu")

# parameters (ticks are set later.)
num_bands = 6
E_Fermi = 22.4023

# x -> the k-path
x = data[:, 0]
# y -> the energy
y = data[:, 1] - E_Fermi

# reshape to [num_bands, num_kpoints]
x = x.reshape([num_bands,-1])
y = y.reshape([num_bands,-1])

# Set up the plot canvas
fig, ax = plt.subplots()

# set xrange to teh first and last k-point
ax.set_xlim([x[0,0],x[0,-1]])

# plot bands
for bands in range(num_bands):
    ax.plot(x[bands], y[bands], linestyle='-', linewidth=1, color="tab:blue")

# setting custom x-ticks
# note that "U|K" adds one additional point.
xticks_positions = [x[0,0], x[0,30], x[0,60], x[0,91], x[0,121],
                    x[0,151],x[0,181]]
xticks_labels = ["Γ", "X", "U|K", "Γ", "L", "W", "X"]
plt.xticks(xticks_positions, xticks_labels)

# Draw vertical grid lines at the high symmetry points
for xc in xticks_positions:
    plt.axvline(x=xc, color='gray', linestyle='--', linewidth=1.0)

# Setting labels and title
plt.ylabel("Energy (eV)")
plt.title("Aluminium Electronic Band Structure")

# Saving the plot
plt.savefig("Al_bands.png", dpi=300)

# Optionally display the plot window
# plt.show()
