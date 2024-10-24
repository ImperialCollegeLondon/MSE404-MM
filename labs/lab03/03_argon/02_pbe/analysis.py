import numpy as np
import matplotlib.pyplot as plt
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
#                               This script checks for converged results. Change the convergence parameter if needed.                        #
#                                                                                                                                            #
#                                               How to run: python3 convergence_processing.py                                                #
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################

def main():
    filename = "data.txt"

    edata = np.loadtxt(filename, delimiter=' ')
    a, etot = edata[:, 0], edata[:, 1]

    minimum_energy = min(etot)
    index_of_minimum = np.where(etot == minimum_energy)
    a_of_minimum_energy = a[index_of_minimum][0]

    print(f"a = {a_of_minimum_energy} gives the minimum energy of {minimum_energy}\n")
    print("You are advised to do a tighter search around this point if you havent already done so")

    plt.figure(figsize=(8, 6))
    plt.scatter(a , etot, color='black', marker='o')

    plt.ylabel("Total Energy (eV)")
    plt.xlabel("Dimer Distance (Å)")
    plt.title("PBE Dimer Distance")
    plt.show()

if __name__ == "__main__":
    main()
