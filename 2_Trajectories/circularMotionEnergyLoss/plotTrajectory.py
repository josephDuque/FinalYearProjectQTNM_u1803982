import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import sys
#from scipy import signal
plt.rcParams.update({'font.size': 13})
plt.rcParams['agg.path.chunksize'] = 100000

pd.options.display.max_rows = 100000000

dataFile = "energyLossCircTrajectoryDeltaT=1.000e-08.tsv"

# Read initial lines - Unstructured data
with open(dataFile) as f:
	reader = csv.reader(f)
	#header = [ next(reader) for x in range(11)] #header array if there are any header rows
f.close()
# Read structured data
powData = pd.read_csv(dataFile, delim_whitespace=True, header=None)
powData.to_numpy()
columns = ["x", "y"]
powData = powData.iloc[:,:-3]
powData.columns = columns[:]

xpos = powData.iloc[100000:100100,0].values
ypos = powData.iloc[100000:100100,1].values

plt.plot(xpos,ypos, linestyle='-', linewidth=0.1, marker=' ', markersize=0.1)
plt.grid()
plt.xlabel("x", fontsize=10)
plt.ylabel("y", fontsize=10)
plt.title("Circular motion with energy loss")
plt.savefig("energyLossCircTrajectory2DeltaT=1.000e-8.png", dpi=300)
plt.show()

"""
ax = powData.plot(x="x",y="y", legend=False, fontsize=9.5, linestyle='none', marker='o', markersize=0.7);
ax.grid()
ax.set_xlabel("x", fontsize=5)
ax.set_ylabel("y", fontsize=5)
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1e'))
plt.tight_layout()

plt.title("Circular motion with energy loss")
plt.savefig("energyLossCircTrajectoryDeltaT=1.000e+4.png", dpi=300)
plt.show()



dataFile = "energyLossCircTrajectoryDeltaT.tsv"

def readFile(dataFile):
        fileObj = open(dataFile, "r") #opens the file in read mode
        words = fileObj.read().splitlines() #puts the file into an array
        fileObj.close()
        return words
"""