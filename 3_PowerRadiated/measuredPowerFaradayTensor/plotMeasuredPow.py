import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
#from scipy import signal
plt.rcParams.update({'font.size': 13})

print(sys.executable)

dataFile = "relMeasuredPowerDeltaT=1.000e-09.tsv"

# Read initial lines - Unstructured data
with open(dataFile) as f:
	reader = csv.reader(f)
	#header = [ next(reader) for x in range(11)] #header array if there are any header rows
f.close()
# Read structured data
powData = pd.read_csv(dataFile, delim_whitespace=True, header=None)
powData.to_numpy()
columns = ["Phase", "Power (W)"]
powData = powData.iloc[:]
powData.columns = columns[:]

ax = powData.plot(x="Phase",y="Power (W)", legend=False, fontsize=12);

ax.grid()
ax.set_ylabel("Power (W)")
plt.tight_layout()\
    
plt.title("dt=1.000e-09")
plt.savefig("relativisticPowerdt=1.000-9.png", dpi=300)
plt.show()