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
    ecut, etot = edata[:, 0], edata[:, 1]

    convergence_parameter = 0.01 #in eV
    print(f"Convergence defined as within {convergence_parameter} eV of the most accurate result")
    flag = 0 # Flag for arrow

    print("ecut (Ry)", " ", "∆_last (meV)")
    print("-----------------------------------------------")
    for row in edata:
        diff = abs(abs(row[1])-abs(etot[-1]))*1000
        if (diff  <= convergence_parameter*1000 and flag==0):
            print(row[0],"  ", diff, "      <-------------")
            flag = 1
        else:
            print(row[0],"  ", diff)

    for i in range(0, len(ecut)):
        if abs(etot[i] - etot[-1]) <= convergence_parameter:
            value = ecut[i]
            print("")
            print(f"Accuracy of {convergence_parameter*1000} meV")
            print(f"Convergence at ecutwfc = {value} Ry")
            break
        else:
            continue

    #plt.figure(figsize=(8, 6))
    #plt.scatter(ecut , etot, color='black', marker='o')

    #plt.ylabel("Total Energy (eV)")
    #plt.xlabel("Energy Cutoff (Ry)")
    #plt.title("Convergence Test")
    #plt.show()

if __name__ == "__main__":
    main()
