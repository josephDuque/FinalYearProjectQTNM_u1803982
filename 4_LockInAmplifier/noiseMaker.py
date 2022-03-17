import numpy as np
from scipy import signal

signalFrequency = 27010796060
signalTimePeriod = 1/signalFrequency
numberOfPeriods = 5000000
samplingRate = 1e11
times = np.arange(0, numberOfPeriods*signalTimePeriod, step=1/samplingRate)

bandwidth = 51e6
T = 4.015
k = 1.38064852e-23

noiseAmplitude = np.sqrt(k*T*bandwidth) #square root to the get voltage amplitude
print(noiseAmplitude)
noise = noiseAmplitude*np.random.normal(size=len(times))

with open("noiseValues.tsv", "w") as txt_file:
    for i in noise:
        txt_file.write(str(i) + "\n")
    
