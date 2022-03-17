import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
#from scipy import signal
plt.rcParams.update({'font.size': 13})

print(sys.executable)

dataFile = "relMeasuredPowerDeltaT=1.000e-10.tsv"

# Read initial lines - Unstructured data
with open(dataFile) as f:
	reader = csv.reader(f)
	#header = [ next(reader) for x in range(11)] #header array if there are any header rows
f.close()
# Read structured data
powData = pd.read_csv(dataFile, delim_whitespace=True, header=None)
powData.to_numpy()
columns = ["Time", "Phase", "Power (W)"]
powData = powData.iloc[:]
powData.columns = columns[:]

#ax = powData.plot(x="Phase",y="Power (W)", legend=False, fontsize=12);


x = powData.iloc[:,0]
y = powData.iloc[:,2]

#ax.grid()
#ax.set_ylabel("Power (W)")

plt.grid()
plt.plot(x,y)
plt.xlabel("Time [s]")
plt.ylabel("Power [W]")
plt.tight_layout()
    

plt.savefig("relativisticPowerdt=1.000e-10.png", dpi=300)
plt.show()