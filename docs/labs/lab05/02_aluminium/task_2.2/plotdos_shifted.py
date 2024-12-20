import numpy as np
import matplotlib.pyplot as plt

# Load the data
data = np.loadtxt('pwscf.dos')

# Adjust the x-values (Energy)
energy = data[:, 0] - 22.402
dos = data[:, 1]
idos = data[:, 2]

# Create the plot
fig, ax1 = plt.subplots()

# Plot Density of States on primary y-axis
ax1.plot(energy, dos, 'b-', label='Density of States')
ax1.set_xlabel('Energy (eV)')
ax1.set_ylabel('Density of states', color='b')
ax1.tick_params(axis='y', labelcolor='b')

# Create secondary y-axis and plot Integrated Density of States
ax2 = ax1.twinx()
ax2.plot(energy, idos, 'r-', label='Integrated Density of States')
ax2.set_ylabel('Integrated density of states', color='r')
ax2.tick_params(axis='y', labelcolor='r')

# set plot range
ax1.set_xlim(energy.min(), energy.max())
ax2.set_xlim(energy.min(), energy.max())

# Set legend
fig.legend(loc='upper left', bbox_to_anchor=(0.12,0.88))

# Format the plot
plt.title('Density of States and Integrated Density of States')

# Save the plot as a PDF
#  plt.savefig('C_diamond_dos.pdf', format='pdf')

# Save the plot as a png
plt.savefig('Al_dos.png', dpi=300)

# Optionally display the plot
# plt.show()
