import numpy as np
import matplotlib.pyplot as plt

# Load data
data_up = np.loadtxt("bands.out.up.gnu")
data_dn = np.loadtxt("bands.out.dn.gnu")

# parameters
num_bands = 16
E_fermi = 12.7648

# Adjust y-values
x_up = data_up[:, 0]
y_up = data_up[:, 1] - E_fermi
# reshape to [num_bands, num_kpoints]
x_up = x_up.reshape([num_bands,-1])
y_up = y_up.reshape([num_bands,-1])

x_dn = data_dn[:, 0]
y_dn = data_dn[:, 1] - E_fermi
# reshape to [num_bands, num_kpoints]
x_dn = x_dn.reshape([num_bands,-1])
y_dn = y_dn.reshape([num_bands,-1])

# Create the plot
fig, ax = plt.subplots()
# set xrange to teh first and last k-point
ax.set_xlim([x_up[0,0],x_up[0,-1]])

# plot bands
for bands in range(num_bands):
    #  ax.plot(x[bands], y[bands], linestyle='-', linewidth=1, color="tab:blue")
    # Plot spin up in red and spin down in blue
    ax.plot(x_up[bands], y_up[bands], color='red', linestyle='-', linewidth=1, label='Spin Up')
    ax.plot(x_dn[bands], y_dn[bands], color='blue', linestyle='-', linewidth=1, label='Spin Down')

# Set x-ticks with custom positions and labels
xticks_positions = [0.0000, 1.732, 2.956, 3.663, 4.163, 5.663, 6.1639]
xticks_labels = ['Γ', 'H', 'N', 'Γ', 'P', 'H|P', 'N']
plt.xticks(xticks_positions, xticks_labels)

# Set the y-axis range
plt.ylim([-10, 15])

# Draw vertical grid lines
for xc in xticks_positions:
    plt.axvline(x=xc, color='gray', linestyle='--', linewidth=1.0)

# Add a horizontal line at the Fermi level
plt.axhline(y=0, color='black', linestyle='--', linewidth=1.0)

# Set labels and title
#  plt.xlabel('Wave Vector')
plt.ylabel('Energy (eV)')
plt.title('Iron Electronic Band Structure')

# Hide the legend (though it's created here for reference)
#ax.legend().set_visible(False)

# Save the plot as a PDF
plt.savefig('Iron_bands.png', dpi=300)
#  plt.savefig('Iron_bands.pdf', format='pdf')

# Optionally display the plot
#plt.show()

