import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import sys
#from scipy import signal
plt.rcParams.update({'font.size': 13})

dataFile = "idealTrajectoryDeltaT=1.000e-09.tsv"

# Read initial lines - Unstructured data
with open(dataFile) as f:
	reader = csv.reader(f)
	#header = [ next(reader) for x in range(11)] #header array if there are any header rows
f.close()
# Read structured data
powData = pd.read_csv(dataFile, delim_whitespace=True, header=None)
powData.to_numpy()
columns = ["phase", "x", "y", "z", "v_x", "v_y", "v_z", "a_x", "a_y", "a_z"]
powData = powData.iloc[:,:]
x = powData.iloc[:,1]
y = powData.iloc[:,2]

v_x = powData.iloc[:,4]
v_y = powData.iloc[:,5]

a_x = powData.iloc[:,7]
a_y = powData.iloc[:,8]
"""
powData.columns = columns[:]
ax = powData.plot(x="x",y="y", legend=False, fontsize=12, linestyle='None', marker='o', markersize=2);

ax.grid()
ax.set_ylabel("y")
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1e'))
"""
plt.plot(x,y)
#plt.plot(v_x, v_y)
#plt.plot(a_x, a_y)

plt.tight_layout()
    
plt.title("Ideal circular motion")
plt.savefig("IdealMotionDt=1.000e-9.png", dpi=300)
plt.show()