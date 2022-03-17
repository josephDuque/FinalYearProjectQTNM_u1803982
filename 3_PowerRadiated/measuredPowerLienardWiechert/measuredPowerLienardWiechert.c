#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#define _USE_MATH_DEFINES
#include<math.h>
#include "hertzianDipole.h"


typedef enum { false, true } bool;

//calculate the fraction of electron speed to light speed from the relativistic kinetic energy:
double calculateBeta(double T) {
  double eMass = 9.10938356*pow(10.0,-31.0);
  double c = 299792458.0;
  
  double factor = pow(T,2.0)/(pow(eMass,2.0)*pow(c,4.0)) + (2.0*T)/(eMass*pow(c,2.0)) + 1;
  double beta = sqrt(1.0 - 1.0/factor);
  return beta;
}

//calculate distance between field point and the moving electron position:
double calculateR(double *fieldPoint, double *electronPos) {
  double distance = sqrt(pow((fieldPoint[0] - electronPos[0]),2.0) + 
		                     pow((fieldPoint[1] - electronPos[1]),2.0) + 
		        	           pow((fieldPoint[2] - electronPos[2]),2.0));
  return distance;
}
//calculate vector between field point and the moving electron position:
void calculateVec(double *fieldPoint, double *electronPos, double *vector) {
  vector[0] = fieldPoint[0] - electronPos[0];
  vector[1] = fieldPoint[1] - electronPos[1];
  vector[2] = fieldPoint[2] - electronPos[2];
}

//calculate unit vector from existing vector
void calculateUnitVec(double *vec, double *unitVec) {
  double magnitude = sqrt(pow(vec[0],2.0) + pow(vec[1],2.0) + pow(vec[2],2.0));
  unitVec[0] = vec[0]/magnitude;
  unitVec[1] = vec[1]/magnitude;
  unitVec[2] = vec[2]/magnitude;
}

//return the correctly scaled vector given its unit vector and magnitude:
void returnProperSizedVec(double *unitVec, double magnitude, double *vec) {
  vec[0] = unitVec[0]*magnitude;
  vec[1] = unitVec[1]*magnitude;
  vec[2] = unitVec[2]*magnitude;
}

//calculate acceleration vector from electron position:
void calculateAccelVec(double *vec, double *accelVec) {
  accelVec[0] = -vec[0];
  accelVec[1] = -vec[1];
  accelVec[2] = -vec[2];
}

//calculate the angle between the acceleration vector (pointing towards the centre) and 
//the unit vector pointing towards the observer.
double calculateTheta(double *accelVec, double *unitVec) {
  double dotProduct = accelVec[0] * unitVec[0] + 
                      accelVec[1] * unitVec[1] +
                      accelVec[2] * unitVec[2] ;
  double accelMag = sqrt(pow(accelVec[0],2.0) + pow(accelVec[1],2.0) + pow(accelVec[2],2.0));
  double theta = acos(dotProduct/accelMag);
  return theta;
}

//magnitude of measured poynting vector for non-relativistic speeds:
double nonRelMeasuredPoyntingVecMag(double R, double angle, double linearV, double effArea, double B) {
  double mu0 = 1.25663706*pow(10.0,-6.0) ;
  double epsilon0 = 8.8541878128*pow(10.0,-12.0);
  double c = 299792458.0;
  double e = 1.602176634*pow(10.0,-19.0);
  double eMass = 9.10938356*pow(10.0,-31.0);
  double radius = (eMass*linearV)/(e*B); //Lamor radius/gyroradius/cyclotron radius
  double centripetalAccel = pow(linearV, 2.0)/radius;

  double factor1 = 1.0/(mu0*c);
  double factor2 = (e/(4.0*M_PI*epsilon0*pow(c,2.0)*R))*centripetalAccel;
  double poyntingMag = factor1*pow(factor2,2.0)*pow(sin(angle),2.0)*effArea;
  return poyntingMag;
}

//magnitude of measured poynting vector for relativistic speeds:
double relMeasuredPoyntingVecMag(double R, double angle, double linearV, double effArea, double B) {
  double mu0 = 1.25663706*pow(10.0,-6.0);
  double epsilon0 = 8.8541878128*pow(10.0,-12.0);
  double c = 299792458.0;
  double e = 1.602176634*pow(10.0,-19.0);
  double eMass = 9.10938356*pow(10.0,-31.0);
  double beta = linearV/c;
  double lorentzFactor = 1.0/pow((1.0 - pow(beta,2.0)),0.5);
  double radius = (lorentzFactor*eMass*linearV)/(e*B); //Lamor radius/gyroradius/cyclotron radius for relativistic motion
  //double centripetalAccel = (pow(lorentzFactor,2.0)*pow(linearV, 2.0))/radius;
  double centripetalAccel = (e*linearV*B)/(lorentzFactor*eMass);
  
  double factor1 = 1.0/(mu0*c);
  double factor2 = (e/(4.0*M_PI*epsilon0*pow(c,2.0)*R))*centripetalAccel;
  double factor3 = 1.0/pow((1.0 - beta*sin(angle)),3.0);
  double factor4 = 1.0 - pow(cos(angle),2.0)/(pow(lorentzFactor,2.0)*pow((1.0 - beta*sin(angle)),2.0));
  double poyntingMag = factor1*pow(factor2,2.0)*factor3*factor4*effArea;
  return poyntingMag;
}

int main() { 
  double c = 299792458.0;
  //double totalTime = 51*5.0e-5;
  double dt = 1.0e-10; //timestep
  int N = 200; //Number of points to calculate
  double B = 1.0; //Idealised magnetic field
  double effArea = 1.0; //effective area of antenna, for idealised purposes set to 1
  //int N = totalTime/dt;
  double e = 1.602176634*pow(10.0,-19.0);
  double eMass = 9.10938356*pow(10.0,-31.0);
  double kineticEnergy = 18.575 * 1000.0 * e; //energy in Joules
  double frequency = (e*B*pow(c,2.0))/(2*M_PI*(eMass*pow(c,2.0) + kineticEnergy)); //frequency, NOT angular frequency
  printf("Frequency calculated: %.4e\n", frequency);
  double wavelength = c/frequency;
  printf("Wavelength: %.4e\n", wavelength);
  //double kineticEnergy =  100.0 * 1000.0 * e; //very relativistic electron
  //double kineticEnergy = 2500 * e; //beta 0.1 threshold is ~ 2578 electronVolts
  double *trajectoryPos, *trajectoryVel, *trajectoryAccel, *observerPos; //Input an electron trajectory, velocity and acceleration
                                       // and initialise observer position
  trajectoryPos = (double *)malloc(3*sizeof(double));
  trajectoryVel = (double *)malloc(3*sizeof(double));
  trajectoryAccel = (double *)malloc(3*sizeof(double));
  double observePoint[3] = {0.0, 5, 0.0};
   
  double *accelVec, *observeVec, *unitObserveVec, *unitAccelVec;
  accelVec = (double *)malloc(3*sizeof(double));
  observeVec = (double *)malloc(3*sizeof(double));
  unitObserveVec = (double *)malloc(3*sizeof(double));
  unitAccelVec = (double *)malloc(3*sizeof(double));

  double *antennaSeparation;
  antennaSeparation = (double *)malloc(3*sizeof(double));
  
  double R; //distance of field point to electron --- will have to vary
  double angle; //angle between acceleration vector and unit vector pointing to observer

  bool determineBeta = true;

  double *phase, *measuredPow;
  double *time;
  time = (double *)malloc(N*sizeof(double));
  phase = (double *)malloc(N*sizeof(double));
  measuredPow = (double *)malloc(N*sizeof(double));

  int i; 
  //determine whether to calculate beta from the relativistic kinetic energy, to choose
  //whether to use either the relativistic or non-relativistic equations:
  if (determineBeta) {
    double beta = calculateBeta(kineticEnergy);
    double electronVel = beta*c;
    double lorentzFactor = 1.0/pow((1.0 - pow(beta,2.0)),0.5);
    //If the electrons velocity is larger than 1/10th the speed of light, use the
    //relativistic calculation, otherwise, use the non-relativistic calculation:
    if (beta > 0.1) {
      printf("Using relativitic calculations...\n");
      printf("velocity calculated: %.4f\n", electronVel);
      double radius = (lorentzFactor*eMass*electronVel)/(e*B); //Lamor radius/gyroradius/cyclotron radius for relativistic motion
      double angularV = electronVel/radius;
      printf("angular velocity: %.4ef\n", angularV);
      double accelMag = (e*electronVel*B)/(eMass);
      printf("lorentzFactor: %.4e\n", lorentzFactor);
      printf("radius: %.4e\n", radius);
      //read trajectory file:
      /*char fileName[40] = "idealTrajectoryDeltaT=";
      char fileType[5] = ".tsv";
      char cutoffVal[15];
      snprintf(cutoffVal, 15, "%.3e", dt);
      char *fileNameExtension;
      fileNameExtension = malloc(strlen(fileName)+1+20);
      strcpy(fileNameExtension, fileName);
      strcat(fileNameExtension, cutoffVal);
      strcat(fileNameExtension, fileType);*/
      char fileNameExtension[40] = "logTime.tsv";
      FILE *posData = fopen(fileNameExtension, "r");
      for (i=0;i<N;i++) {
        //fscanf(posData, "%lf \t %lf \t %lf \n", &trajectoryPos[0], &trajectoryPos[1], &trajectoryPos[2]);
        //phase[i] = angularV*dt*i;
        fscanf(posData, "%lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf f\n",
                       &time[i], //&phase[i],
                       &trajectoryPos[0], &trajectoryPos[1], &trajectoryPos[2],
                       &trajectoryVel[0], &trajectoryVel[1], &trajectoryVel[2]);
                       //&trajectoryAccel[0], &trajectoryAccel[1], &trajectoryAccel[2]);
        //printf("position: %.3e \t %.3e \t %.3e\n", trajectoryPos[0], trajectoryPos[1], trajectoryPos[2]);
        R = calculateR(observePoint, trajectoryPos);
        //printf("distance calculated: %.3e\n", R);
        calculateVec(observePoint, trajectoryPos, observeVec);  //finds vector pointing to observer
        calculateVec(trajectoryPos, observePoint, antennaSeparation);  //finds vector pointing to observer
        //printf("observing vector: %.3e \t %.3e \t %.3e\n", observeVec[0], observeVec[1], observeVec[2]);
        calculateUnitVec(observeVec, unitObserveVec);  //finds unit vector pointing to observer
        //printf("observing unit vector: %.3e \t %.3e \t %.3e\n", unitObserveVec[0], unitObserveVec[1], unitObserveVec[2]);
        calculateAccelVec(trajectoryPos, accelVec);  //finds acceleration vector direction
        calculateUnitVec(accelVec, unitAccelVec);   //normalises direction
        returnProperSizedVec(unitAccelVec, accelMag, accelVec); //scales acceleration vector correctly
        //printf("acceleration vector: %.3e \t %.3e \t %.3e\n", accelVec[0], accelVec[1], accelVec[2]);
        effArea = calcEffArea(antennaSeparation, wavelength);
        //printf("Effective area: %.4e\n", effArea);
        angle = calculateTheta(accelVec, unitObserveVec);  //finds angle between acceleration vector and unit vector to observer
        //printf("Angle calculated: %.3e\n", angle);
        measuredPow[i] = relMeasuredPoyntingVecMag(R,angle,electronVel,effArea,B);
      }
      //free(fileNameExtension);
      fclose(posData);

      char fileName1[40] = "relMeasuredPowerDeltaT=";
      char fileType1[5] = ".tsv";
      char cutoffVal1[15];
      snprintf(cutoffVal1, 15, "%.3e", dt);
      char *fileNameExtension1;
      fileNameExtension1 = malloc(strlen(fileName1)+1+20);
      strcpy(fileNameExtension1, fileName1);
      strcat(fileNameExtension1, cutoffVal1);
      strcat(fileNameExtension1, fileType1);
  
      FILE *powData = fopen(fileNameExtension1, "w");
      for (i=0;i<N;i++) {
        fprintf(powData, "%.3e \t %.3e \t %.3e \n", time[i], phase[i], measuredPow[i]*80000);
      }
      fclose(powData);
      free(fileNameExtension1);
    }
    else {
      printf("Using non-relativitic calculations...");
      printf("velocity calculated: %.4f\n", electronVel);
      double electronVelTheoretical = sqrt((2.0*kineticEnergy)/eMass);
      double radius = (eMass*electronVelTheoretical)/(e*B); //Lamor radius/gyroradius/cyclotron radius
      double angularV = electronVelTheoretical/radius;
      
      char fileName[40] = "idealTrajectoryDeltaT=";
      char fileType[5] = ".tsv";
      char cutoffVal[15];
      snprintf(cutoffVal, 15, "%.3e", dt);
      char *fileNameExtension;
      fileNameExtension = malloc(strlen(fileName)+1+20);
      strcpy(fileNameExtension, fileName);
      strcat(fileNameExtension, cutoffVal);
      strcat(fileNameExtension, fileType);
      FILE *posData = fopen(fileNameExtension, "r");
      for (i=0;i<N;i++) {
        fscanf(posData, "%lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \n",
                       &time[i], &phase[i],
                       &trajectoryPos[0], &trajectoryPos[1], &trajectoryPos[2],
                       &trajectoryVel[0], &trajectoryVel[1], &trajectoryVel[2],
                       &trajectoryAccel[0], &trajectoryAccel[1], &trajectoryAccel[2]);
        //printf("position: %.3e \t %.3e \t %.3e\n", trajectoryPos[0], trajectoryPos[1], trajectoryPos[2]);
        R = calculateR(observePoint, trajectoryPos);
        //printf("distance calculated: %.3e\n", R);
        calculateVec(observePoint, trajectoryPos, observeVec);  //finds vector pointing to observer
        calculateVec(observePoint, trajectoryPos, observeVec);  //finds vector pointing to observer
        //printf("observing vector: %.3e \t %.3e \t %.3e\n", observeVec[0], observeVec[1], observeVec[2]);
        calculateUnitVec(observeVec, unitObserveVec);  //finds unit vector pointing to observer
        //printf("observing unit vector: %.3e \t %.3e \t %.3e\n", unitObserveVec[0], unitObserveVec[1], unitObserveVec[2]);
        calculateAccelVec(trajectoryPos, accelVec);  //finds acceleration vector
        //printf("acceleration vector: %.3e \t %.3e \t %.3e\n", accelVec[0], accelVec[1], accelVec[2]);
        angle = calculateTheta(accelVec, unitObserveVec);  //finds angle between acceleration vector and unit vector to observer
        //printf("Angle calculated: %.3e\n", angle);
        measuredPow[i] = nonRelMeasuredPoyntingVecMag(R,angle,electronVel,effArea,B);
      }
      free(fileNameExtension);
      fclose(posData);

      char fileName1[40] = "nonRelMeasuredPowerDeltaT=";
      char fileType1[5] = ".tsv";
      char cutoffVal1[15];
      snprintf(cutoffVal1, 15, "%.3e", dt);
      char *fileNameExtension1;
      fileNameExtension1 = malloc(strlen(fileName1)+1+20);
      strcpy(fileNameExtension1, fileName1);
      strcat(fileNameExtension1, cutoffVal1);
      strcat(fileNameExtension1, fileType1);
      FILE *powData = fopen(fileNameExtension1, "w");
      for (i=0;i<N;i++) {
        fprintf(powData, "%.3e \t %.3e \n", phase[i], measuredPow[i]);
      }
      fclose(powData);
      free(fileNameExtension1);
    }
  }
  //use non-relativisitc equation for any energy:
  else {
    double electronVel = sqrt((2.0*kineticEnergy)/eMass);
    double radius = (eMass*electronVel)/(e*B); //Lamor radius/gyroradius/cyclotron radius
    double angularV = electronVel/radius;
    printf("velocity calculated: %.4f\n", electronVel);
    for (i=0;i<N;i++) {
      time[i] = i*dt;
      phase[i] = angularV * time[i];
      measuredPow[i] = nonRelMeasuredPoyntingVecMag(R,time[i],electronVel,effArea,B);
    }
    char fileName[40] = "nonRelMeasuredPowerDeltaT=";
    char fileType[5] = ".tsv";
    char cutoffVal[15];
    snprintf(cutoffVal, 15, "%.1e", dt);
    char *fileNameExtension;
    fileNameExtension = malloc(strlen(fileName)+1+20);
    strcpy(fileNameExtension, fileName);
    strcat(fileNameExtension, cutoffVal);
    strcat(fileNameExtension, fileType);
  
    FILE *powData = fopen(fileNameExtension, "w");
    for (i=0;i<N;i++) {
      fprintf(powData, "%.3e \t %.3e \n", phase[i], measuredPow[i]);
    }
    fclose(powData);
    free(fileNameExtension);
  }
  
  free(accelVec); free(observeVec); free(unitObserveVec); free(unitAccelVec); free(antennaSeparation);
  free(time); free(phase); free(measuredPow); free(trajectoryPos); free(trajectoryVel); free(trajectoryAccel);
}