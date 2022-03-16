#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#define _USE_MATH_DEFINES
#include<math.h>

//calculate the fraction of electron speed to light speed from the relativistic kinetic energy:
double calculateBeta(double T) {
  double eMass = 9.10938356*pow(10.0,-31.0);
  double c = 299792458.0;
  
  double factor = pow(T,2.0)/(pow(eMass,2.0)*pow(c,4.0)) + (2.0*T)/(eMass*pow(c,2.0)) + 1;
  double beta = sqrt(1.0 - 1.0/factor);
  return beta;
}

//calculate the speed of the electron from the relativistic kinetic energy:
double calculateV(double T) {
  double eMass = 9.10938356*pow(10.0,-31.0);
  double c = 299792458.0;
  
  double factor = pow(T,2.0)/(pow(eMass,2.0)*pow(c,4.0)) + (2.0*T)/(eMass*pow(c,2.0)) + 1;
  double beta = sqrt(1.0 - 1.0/factor);
  double electronVel = beta*c;
  return electronVel;
}

//finding the phase shift in the equation of motion from the starting position:
double findPhaseShift(double *startingPos) {
  double phase = atan(startingPos[1]/startingPos[0]);
  return phase;
}

//find electron position from circular motion equation:
void findPos(double *electronPos, double radius, double phase) {
  electronPos[0] = radius*cos(phase);
  electronPos[1] = radius*sin(phase);
  electronPos[2] = 0;
}
//find electron velocity
void findVel(double *electronVel, double radius, double phase, double angVel) {
  electronVel[0] = -angVel*radius*sin(phase);
  electronVel[1] = angVel*radius*cos(phase);
  electronVel[2] = 0;
}
//find electron acceleration
void findAccel(double *electronAccel, double radius, double phase, double angVel) {
  electronAccel[0] = -pow(angVel,2.0)*radius*cos(phase);
  electronAccel[1] = -pow(angVel,2.0)*radius*sin(phase);
  electronAccel[2] = 0;
}


int main() {
  double c = 299792458.0;
  double B = 1.0;
  double e = 1.602176634*pow(10.0,-19.0);
  double eMass = 9.10938356*pow(10.0,-31.0);
  double kineticEnergy = 18.575 * 1000.0 * e; //kinetic energy in Joules
  double beta = calculateBeta(kineticEnergy);
  double lorentzFactor = 1.0/pow((1.0 - pow(beta,2.0)),0.5);
  double linearV = beta*c;
  double radius = (lorentzFactor*eMass*linearV)/(e*B); //Lamor radius/gyroradius/cyclotron radius for relativistic motion
  
  double dt = 1.66e-11;
  double startingPos[3] = {radius*cos(M_PI/4.0), radius*sin(M_PI/4.0), 0.0};
  int N = 100000000; //Number of points
  //double centripetalAccel = (pow(lorentzFactor,2.0)*pow(linearV, 2.0))/radius;
  double angularV = linearV/radius; //may include lorentz factor
  printf("Linear velocity: %.4e\n", linearV);
  printf("Angular velocity: %.4e\n", angularV);
  printf("Radius: %.4e\n", radius);
  double phase = findPhaseShift(startingPos); 
  double time, realPhase;
  
  double *electronPos, *electronVel, *electronAccel;
  electronPos = (double *)malloc(3*sizeof(double));
  electronVel = (double *)malloc(3*sizeof(double));
  electronAccel = (double*)malloc(3*sizeof(double));
  
  int i,j;
  char fileName[40] = "idealTrajectoryDeltaT=";
  char fileType[5] = ".tsv";
  char cutoffVal[15];
  snprintf(cutoffVal, 15, "%.3e", dt);
  char *fileNameExtension;
  fileNameExtension = malloc(strlen(fileName)+1+20);
  strcpy(fileNameExtension, fileName);
  strcat(fileNameExtension, cutoffVal);
  strcat(fileNameExtension, fileType);

  FILE *posData = fopen(fileNameExtension, "w");
  for (i=0;i<N;i++) {
    time = i*dt;
    realPhase = angularV*time + phase;
    findPos(electronPos, radius, realPhase);
    findVel(electronVel, radius, realPhase, angularV);
    findAccel(electronAccel, radius, realPhase, angularV);
    fprintf(posData, "%.3e \t %.3e \t %.3e \t %.3e \t %.3e \t %.3e \t %.3e \t %.3e \t %.3e \t %.3e \t %.3e \n",
            time, realPhase, 
            electronPos[0], electronPos[1], electronPos[2],
            electronVel[0], electronVel[1], electronVel[2], 
            electronAccel[0], electronAccel[1], electronAccel[2]);
  }
  fclose(posData);
  printf("File saved: %s", fileNameExtension);
  free(fileNameExtension);

  free(electronPos); free(electronVel); free(electronAccel);
}