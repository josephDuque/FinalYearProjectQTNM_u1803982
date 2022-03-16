import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
#from scipy import signal
plt.rcParams.update({'font.size': 13})

dataFile = "bathtubZaxis_solenoidCoil.tsv"

# Read initial lines - Unstructured data
with open(dataFile) as f:
	reader = csv.reader(f)
	#header = [ next(reader) for x in range(11)] #header array if there are any header rows
f.close()
# Read structured data
powData = pd.read_csv(dataFile, delim_whitespace=True, header=None)
powData.to_numpy()
columns = ["z [m]", "Total magnetic field [T]"]
powData = powData.iloc[:]
powData.columns = columns[:]

ax = powData.plot(x="z [m]",y="Total magnetic field [T]", legend=False, fontsize=12);

ax.grid()
ax.set_ylabel("Total magnetic field [T]")
plt.tight_layout()
    
plt.title("Magnetic bathtub trap - 2 solenoids of 11 coils")
plt.savefig("bathtubTrapCoils=11_Width=0,2", dpi=300, bbox_inches="tight")
plt.show()