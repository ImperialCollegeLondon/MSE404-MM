import numpy as np
import matplotlib.pyplot as plt

# Reading the data from the file
data = np.loadtxt("bands.out.gnu")

# Adjust the y-values (second column by subtracting 13.993)
x = data[:, 0]
y = data[:, 1] - 13.993

# Set up the plot
fig, ax = plt.subplots()

ax.plot(x, y, linestyle='-', linewidth=1)

# Setting custom x-ticks
xticks_positions = [0.0000, 1.0607, 1.4142, 2.4142, 3.2802, 4.1463, 4.6463, 5.3534]
xticks_labels = ["Γ", "K", "X", "Γ", "L", "X", "W", "L"]
plt.xticks(xticks_positions, xticks_labels)

# Draw vertical grid lines
for xc in xticks_positions:
    plt.axvline(x=xc, color='gray', linestyle='--', linewidth=1.0)

# Setting labels and title
plt.ylabel("Energy (eV)")
plt.title("Carbon Diamond Electronic Band Structure")

# Hide the plot legend
ax.legend().set_visible(False)

# Saving the plot as a PDF
plt.savefig("C_diamond_bands.pdf", format='pdf')

# Optionally display the plot window
# plt.show()
