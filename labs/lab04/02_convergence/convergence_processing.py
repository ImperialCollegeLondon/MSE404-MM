import numpy as np
import matplotlib.pyplot as plt
 
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#                       This script works as is. An arrow will point towards your converged result.                                   #
#                       To run this code all you need to run is : python3 kpt_convergence_processing.py                               #
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
 
def main():
    filename = "data_kpt.txt"
 
    edata = np.loadtxt(filename, delimiter=' ')
    kpt, etot = edata[:, 0], edata[:, 1]
 
    convergence_parameter = 0.01 #in eV
    flag = 0 # Flag for arrow
 
    print("K Points nxnxn", "       ", "∆_last (meV/atom)")
    print("-----------------------------------------------")
    for row in edata:
        diff = abs(abs(row[1])-abs(etot[-1]))*1000
        if (diff  <= convergence_parameter*1000 and flag==0):
            print(row[0],"                ", diff, "      <-------------")
            flag = 1
        else:
            print(row[0],"                ", diff)
 
    for i in range(0, len(kpt)):
        if abs(etot[i] - etot[-1]) <= convergence_parameter:
            value = kpt[i]
            print("")
            print(f"Accuracy of {convergence_parameter*1000} meV")
            print(f"Convergence at Kpts = {value}")
            break
        else:
            continue
 
    plt.figure(figsize=(8, 6))
    plt.scatter(kpt , etot, color='black', marker='o')
 
    for row in edata:
        diff = abs(abs(row[1])-abs(etot[-1]))*1000
        if (diff  <= convergence_parameter*1000):
            plt.scatter(row[0] , row[1], color='red', marker='o')
 
    plt.ylabel("Total Energy (eV)")
    plt.xlabel("K Points nxnxn")
    plt.ticklabel_format(useOffset=False)
    plt.title("Convergence Testing")
    plt.show()
 
if __name__ == "__main__":
    main()
