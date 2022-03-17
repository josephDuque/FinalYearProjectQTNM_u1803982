#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#define _USE_MATH_DEFINES
#include<math.h>

//ButterLowPass - filters data using a Butterworth low pass filter
double applyFilter(double cutoff, int order, double samplingRate, double waveformDatapoint) {
  double nyquistFreq = 0.5 * samplingRate;
  double normalCutoff = cutoff / nyquistFreq;
  double fraction = waveformDatapoint/normalCutoff;
  double filteredData = 1.0 / sqrt(1.0 + pow(fraction, 2.0*order));
  return filteredData;
}

//lockInAmplifier - caluclates quadrature and prints results to text file.
void lockInAmp(double *signal, double *reference, double *referencePhaseShift,
              double cutoff, int order, double samplingRate, double *magnitudes, double *phases) { 
  int M, i;                       
  M = (int)sizeof(signal)/sizeof(signal[0]);
  printf("M = %d\n", M);
  
  double *mixSigAndRef, *mixSigAndRefShift, *filterMixSig, *filterMixSigShift;

  mixSigAndRef = (double *)malloc(M * sizeof(double));
  mixSigAndRefShift = (double *)malloc(M * sizeof(double));
  filterMixSig = (double *)malloc(M * sizeof(double));
  filterMixSigShift = (double *)malloc(M * sizeof(double));

  for (i=0;i<M;i++) {
    mixSigAndRef[i] = signal[i]*reference[i];
    mixSigAndRefShift[i] = signal[i]*referencePhaseShift[i];
    filterMixSig[i] = applyFilter(cutoff, order, samplingRate, mixSigAndRef[i]);
    filterMixSigShift[i] = applyFilter(cutoff, order, samplingRate, mixSigAndRefShift[i]);
    magnitudes[i] = sqrt(pow(filterMixSig[i],2) + pow(filterMixSigShift[i],2));
    phases[i] = atan2(filterMixSigShift[i], filterMixSig[i]);
  }
 
  free(mixSigAndRef); free(mixSigAndRefShift);
  free(filterMixSig); free(filterMixSigShift);
}

//Function to randomise signal between an interval
double drand ( double low, double high )
{
    return ( (double)rand() * ( high - low ) ) / (double)RAND_MAX + low;
}

int main() {
  double signalMin = 27010796060;
  double signalMax = 27010796060 + 51.0e6;
  double signalAmplitude = 3.16e-8; //estimate to the amplitude of the voltage from a power signal of 1e-15
  double signalFreq = drand(signalMin, signalMax); //signal frequency randomised inside the bandwidth of interest (dont know what it is)
  double amplifierScanBandwidth = 51.0e4; //each amplifier scans a bandwidth of 0.51 MHz of the total bandwidth, so 100 lock-ins

  double referenceFrequency = 27010796060; //reference frequency of an endpoint electron with a zero mass neutrino, this will be shifted for each lock in amplifier
  double referencePhaseShift = (1.0/4.0)*M_PI; //Some phase shift
  double referenceAmplitude = 3.16e-8;  //estimate to the amplitude of the voltage from a power signal of 1e-15

  double signalPeriod = 1.0/referenceFrequency;
  int numPeriods = 5000000;
  double samplingRate = 1.1e11;  //This affects the filter cutoff, should be 4 times larger than signal frequency

  double *times;
  double *signalSin;

  int N, i, j;
  double step = 1.0/samplingRate; //Time step dt
  double totalSignalTime = numPeriods*signalPeriod; 
  N = (totalSignalTime)/step + 1;
  printf("Total signal time = %f\n", totalSignalTime);
  printf("N = %d\n", N);
  times = (double *)malloc(N * sizeof(double));
  signalSin = (double *)malloc(N * sizeof(double));
  //Voltage sine wave of signal 
  for (i=0;i<N;i++) {
    times[i] = i*step;
    signalSin[i] = sin(2.0*M_PI*signalFreq*times[i]) * signalAmplitude;
  }

//Add noise:
  double *noise; double *noisySignal;
  //double mean = times[(int)N/2];
  noise = (double *)malloc(N * sizeof(double));
  //read noise off a tsv file:
  noisySignal = (double *)malloc(N * sizeof(double));
  FILE *noiseFile = fopen("noiseValues.tsv", "r");
  for (i=0;i<N;i++) {
    fscanf(noiseFile, "%lf\n", &noise[i]);
    noisySignal[i] = signalSin[i] + noise[i];
  }
  fclose(noiseFile);


  //Produce reference and 90 degree phase shifted reference
  double *referenceSin, *referenceCos;
  double referenceShift = amplifierScanBandwidth/2.0; //reference frequency at half the bandwidth scanned
  double refFreqUpdated = referenceFrequency + referenceShift;
  referenceSin = (double *)malloc(N * sizeof(double));
  referenceCos = (double *)malloc(N * sizeof(double));
  for (i=0;i<N;i++) {
    referenceSin[i] = sin(2.0*M_PI*refFreqUpdated*times[i] + referencePhaseShift) * referenceAmplitude;
    referenceCos[i] = cos(2.0*M_PI*refFreqUpdated*times[i] + referencePhaseShift) * referenceAmplitude;
  }

//Setting up cutoff frequency of filter
  //double cutoff[3] = {1e5, 1e4, 1e3};
  double cutoff = 1e5;

//Apply lock in apmplifiers using different reference frequencies, and print results to data files:
  double *magnitudes, *phases;

  magnitudes = (double *)malloc(N * sizeof(double));
  phases = (double *)malloc(N * sizeof(double));

  int noOfLockInProbes = (int)51.0e6/51.0e4 ;
  for (i=0;i<noOfLockInProbes;i++) {
    printf("Running through loop: %d\n\n\n", i);
    char fileName[40] = "lockInAmp";
    char fileType[5] = ".tsv";
    char amplifierNo[15];
    snprintf(amplifierNo, 15, "%d", (i+1));
    char *fileNameExtension;
    fileNameExtension = malloc(strlen(fileName)+1+20);
    strcpy(fileNameExtension, fileName);
    strcat(fileNameExtension, amplifierNo);
    strcat(fileNameExtension, fileType);
    
    //  signal with noise, in phase reference, quadrature, cutoff, order, sampling rate, and the outputs
    lockInAmp(noisySignal, referenceSin, referenceCos, cutoff, 2, samplingRate, magnitudes, phases);

    FILE *file = fopen(fileNameExtension, "w");
    for (j=0;j<N;j++) {
      printf("%f \t %.2e \t %.2e \t %.2e \t %.2e\n", times[j], signalSin[j], noisySignal[j], magnitudes[j], phases[j]);
      //fprintf(file, "%f \n", times[j]);
    }
    fclose(file); //\t %f \t %f , magnitudes[j], phases[j]
    refFreqUpdated = referenceFrequency + (i+1)*amplifierScanBandwidth + referenceShift; //change reference frequency for the next lock-in amplifier
    for (j=0;j<N;j++) {
      referenceSin[j] = sin(2.0*M_PI*refFreqUpdated*times[j] + referencePhaseShift) * referenceAmplitude;
      referenceCos[j] = cos(2.0*M_PI*refFreqUpdated*times[j] + referencePhaseShift) * referenceAmplitude;
    }
 
    
    free(fileNameExtension);
  }
  
  free(magnitudes); free(phases);
  free(referenceSin); free(referenceCos);
  free(times); free(signalSin);
  free(noise); free(noisySignal);
  return 0; 
}