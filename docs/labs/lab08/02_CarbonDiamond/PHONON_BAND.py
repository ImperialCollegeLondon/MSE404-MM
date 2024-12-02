import numpy as np
import matplotlib.pyplot as plt

#Labels for the high-symmetry q-points. We have prepared them for you.
qmarker_lab = [r'$\Gamma$',r'$K$',r'$X$',r'$\Gamma$',r'$L$',r'$X$',r'$W$',r'$L$']
#Name of the file containing the coordinates of the q-points and frequencies
freq_file = 'CD-bands.freq.gp'
#Input file for using matdyn.x
input_file_for_matdyn = '04_CD_matdyn-bands.in'

#Parameters for the plot. These are recommended but you are encouraged to play around with them
fontsize = 12.5
figsize = (6.5,5)

###################################### The data files will be read below
f = open(freq_file,'r+')  #Reads the frequency data file
lines = f.readlines()
f.close()

g = open(input_file_for_matdyn,'r+')  #Reads the input file
inputlines = g.readlines()
g.close()

bands = []  #Converting the file into a matrix
for j in range(len(lines)):
    line = lines[j].split()[:]
    bands.append(np.array(line,dtype=np.float64))
bands = np.array(bands)

count = 0   #Reading the coordinates of the high-symmetry points into a list
L_seg = []
for j in range(len(inputlines)):
    #The if conditions below are designed to locate the lines containing the high-symmetry points
    if count == 1:
        count += 1 #This is to skip the line in the input file that contains the number of segments. See below.
    elif count == 2:
        L_seg.append(eval(inputlines[j].split()[-2])) #This starts actually reading the number of points along each segment
    if '/' in inputlines[j]: #Tells the script to look for the '/' sign, the coordinates of the high-symmetry points will be found two lines below this line.
        count += 1

marker_pos = []
for j in range(len(L_seg)):
    marker_pos.append(sum(L_seg[:j]))  #Determine the location of the labels based on the number of q-points between each high-symmetry point
marker_pos = np.array(marker_pos)

qpoints = bands[:,0]
qmarker_pos = qpoints[marker_pos]

#Start making the plots below
plt.figure(figsize=figsize)
plt.vlines(qmarker_pos,np.min(bands[:,1:]),np.max(bands[:,1:])*1.05,color='black',linestyle='--',alpha=0.4)
#Read the bands matrix correctly to produce a band plot one by one
for j in range(1,bands.shape[-1]):
    plt.plot(qpoints,bands[:,j],c='purple')
plt.xticks(qmarker_pos,qmarker_lab,fontsize=fontsize)
plt.yticks(fontsize=fontsize)
plt.ylabel(r'$\omega (\mathrm{cm}^{-1})$',fontsize=fontsize)
plt.savefig('bands', dpi=300, bbox_inches="tight",pad_inches=0.01)
plt.show()
