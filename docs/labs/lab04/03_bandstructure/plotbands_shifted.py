import numpy as np
import matplotlib.pyplot as plt

# Reading the data from the file
data = np.loadtxt("bands.out.gnu")

# Adjust the y-values (second column by subtracting 13.993)
x = data[:, 0]
y = data[:, 1] - 13.993
# reshape, we have a total of 8 bands, each line contain 182 points
x = x.reshape([8,182])
y = y.reshape([8,182])

# Set up the plot
fig, ax = plt.subplots()

# set xrange
ax.set_xlim([x[0,0],x[0,-2]])

# plot bands
for bands in range(8):
    ax.plot(x[bands], y[bands], linestyle='-', linewidth=1, color="tab:blue")

#Setting custom x-ticks
xticks_positions = [x[0,0], x[0,31], x[0,61], x[0,91], x[0,121],
                    x[0,151],x[0,181]]
xticks_labels = ["Γ", "X", "U|K", "Γ", "L", "W", "X"]
plt.xticks(xticks_positions, xticks_labels)

# Draw vertical grid lines
for xc in xticks_positions:
    plt.axvline(x=xc, color='gray', linestyle='--', linewidth=1.0)

# Setting labels and title
plt.ylabel("Energy (eV)")
plt.title("Carbon Diamond Electronic Band Structure")

# Hide the plot legend
#  ax.legend().set_visible(False)

# Saving the plot
#  plt.savefig("C_diamond_bands.pdf", format='pdf')
plt.savefig("C_diamond_bands.png", dpi=300)

# Optionally display the plot window
# plt.show()
