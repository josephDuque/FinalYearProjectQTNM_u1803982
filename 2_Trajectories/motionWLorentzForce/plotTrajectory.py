#!/usr/local/bin/python3
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 13 05:32:46 2022

@author: josephDuque
"""
import sys, getopt
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.cbook as cbook
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
from scipy import ndimage
# ----------------------------------------------------------
plt.rcParams.update({
    "text.usetex": True,
    "font.size"  : 12,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica"]})
# for Palatino and other serif fonts use:
plt.rcParams.update({
    "text.usetex": True,
    "font.size"  : 12,
    "font.family": "serif",
    "font.serif": ["Palatino"],
})
# ----------------------------------------------------------
#print("backend", plt.rcParams["backend"])
#plt.rcParams["backend"] = "TkAgg" # doesn't actually set the backend
#matplotlib.use("TkAgg")
print("backend", plt.rcParams["backend"])
#print("sys.float_info.dig = ", sys.float_info.dig)
#print("sys.float_info.mant_dig = ", sys.float_info.mant_dig)
lowEps = np.finfo(float).eps*np.float(100.0)
#print("lower precision limit=",lowEps)
# ----------------------------------------------------------

def main(argv):
    inputfile = ''
    outputfile = ''
    try:
       opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
    except getopt.GetoptError:
       print('plotTrajectory.py -i <dataToPlotDile> -o <hardCopyFileName>');
       sys.exit(2)
    for opt, arg in opts:
       if opt == '-h':
          print('plotTrajectory.py -i <dataToPlotDile> -o <hardCopyFileName>');
          sys.exit()
       elif opt in ("-i", "--ifile"):
          inputfile = arg
       elif opt in ("-o", "--ofile"):
          outputfile = arg
    print('File with data to plot <Input file> is "', inputfile)
    print('Filename of hard copy file <Output file> is "', outputfile)
    # --- Reading data 
    t0,posx,posy,posz,ux,uy,uz = \
        np.loadtxt(open(inputfile,'rt').readlines()[:-1], delimiter='\t', skiprows=12, unpack=True);
    # --- Creating plot
    #fig = plt.figure()
    #ax = fig.add_subplot(projection='3d')
    #scatter = ax.scatter(posx,posy,posz,c=posz,cmap='viridis',alpha=0.75)
    # legend1 = ax.legend(*scatter.legend_elements(),
    #                     loc="upper left", title=r"$z-$Position of a Cyclotron Trajectory",fontsize=12)
    #ax.add_artist(legend1)
    
    fig = plt.figure(figsize=(6,4.5), dpi=144)
    ax = Axes3D(fig);
    line = plt.plot(posx,posy,posz,lw=0.2,c='k')[0]
    ax.view_init(azim=44.5,elev=15.)

    ax.grid(True)
    ax.set_xlabel(r'$x-$Position',fontsize=12)
    ax.set_ylabel(r'$y-$Position',fontsize=12);
    ax.set_zlabel(r'$z-$Position',fontsize=12);
    
    # string01 = '$\max(\phi) = ' + str(np.amax(col07)) + '$' 
    # string02 = '$\min(\phi) = ' + str(np.amin(col07)) + '$' 
    # ax.text(2, 80.0, r"$\Pi_1 = \frac{p\,\rho\,c_0^2}{\mu}$", color="k", fontsize=18)
    # ax.text(2, 74.0, r"$\Pi_2 = \frac{p\,c_0}{h_0}$", color="k", fontsize=18)
    # ax.text(2, 68.0, r"$\mu = 1.3e-2, \, \rho=1005.0$", color="k", fontsize=18)
    # ax.text(2, 60.0, string01, color="k", fontsize=18)
    # ax.text(2, 55.0, string02, color="k", fontsize=18)
    
    plt.show()


# -----------------------------------------------------------------
# Main function call
if __name__ == "__main__":
    if (len(sys.argv)>1):
        main(sys.argv[1:]);
    else:
        print('Please provide input file to plot.')
        print('plotTrajectory.py -i <dataToPlotDile> -o <hardCopyFileName>');

        
# End of main function call
# -----------------------------------------------------------------
