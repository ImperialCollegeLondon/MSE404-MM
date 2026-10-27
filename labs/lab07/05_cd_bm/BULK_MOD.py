#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#      	 This is a plotting code. Put your data into data.txt in column format. <column 1 = %change in volume> <column 2 = total energy>               #
#      	 					To use this code, issue the command:							    #
#      	 						python3 plot.py									    #
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
import numpy as np
import matplotlib.pyplot as plt


def main():
    data = np.loadtxt("data.txt") #Load the data

    vs, etot = data[:,0], data[:,1]

    b2a = 0.529177249 #Ratio of Bohrs to Å
    vs *= b2a**3 #Convert the volume from Bohr^3 to Å^3

    vgrid = np.linspace(min(vs),max(vs),100)

    a, b, c = np.polyfit(vs,etot,deg=2) #numpy.polyfit is a function for fitting a parabola to the data points and extracting the coefficients.

    bm = -2*a*vs[int(len(vs)/2)] #The value of volume at the center of the array is the relaxed unit-cell volume.

    ry2j = 2.18e-18 #Convert Rydbergs to Joules
    ang2m = 1e-10 #Convert Å to meters
    convert_ratio = ry2j/ang2m**3 #Convert the bulk modulus into SI units
    print(f'The bulk modulus is {bm*convert_ratio/10**9:.3e} GPa.')

    fontsize = 15 #Fontsize of the plot
    plt.figure(figsize=(10.5, 6))
    plt.plot(vgrid, c+b*vgrid+a*vgrid**2,c='orange')
    plt.scatter(vs, etot)
    plt.xlabel(r"$V  (\mathrm{\AA}^3)$",fontsize=fontsize)
    plt.ylabel("Total Energy (Ry)",fontsize=fontsize)
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    plt.title("PES against V",fontsize=fontsize)
    plt.savefig('bm.png', dpi=300, bbox_inches="tight",pad_inches=0.01)
    plt.show()

if __name__ == "__main__":
	main()
