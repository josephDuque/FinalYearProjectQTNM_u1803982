#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#define _use_math_defines
#include<math.h>
#define c 299792458.0
#define mu 1.25663706e-6 //permeability of free space
#define e 1.602176634*pow(10.0,-19.0)
#define eMass 9.10938356*pow(10.0,-31.0)

//calculate the fraction of electron speed to light speed from the relativistic kinetic energy:
double calculateBeta(double T) {
  //double eMass = 9.10938356*pow(10.0,-31.0);
  //double c = 299792458.0;
  
  double factor = pow(T,2.0)/(pow(eMass,2.0)*pow(c,4.0)) + (2.0*T)/(eMass*pow(c,2.0)) + 1;
  double beta = sqrt(1.0 - 1.0/factor);
  return beta;
}

//dot product of two 3-vectors:
double dotProduct3vectors(double *vec1, double *vec2){
  int i;
  double product = 0.0;
  for (i=0;i<3;i++){
    product = vec1[i]*vec2[i] + product;
  }
  return product;
}

//dot product of two 4-vectors:
double dotProduct4vectors(double *vec1, double *vec2){
  int i;
  double product = vec1[0]*vec2[0] - (vec1[1]*vec2[1] + vec1[2]*vec2[2] + vec1[3]*vec2[3]);
  return product;
}

//calculate distance between field point and the moving electron position:
double calculateR(double *fieldPoint, double *electronPos) {
  double distance = sqrt(pow((fieldPoint[0] - electronPos[0]),2.0) + 
		                     pow((fieldPoint[1] - electronPos[1]),2.0) + 
		        	           pow((fieldPoint[2] - electronPos[2]),2.0));
  return distance;
}

//calculate unit vector from existing vector
void calculateUnitVec(double *vec, double *unitVec) {
  double magnitude = sqrt(pow(vec[0],2.0) + pow(vec[1],2.0) + pow(vec[2],2.0));
  unitVec[0] = vec[0]/magnitude;
  unitVec[1] = vec[1]/magnitude;
  unitVec[2] = vec[2]/magnitude;
}

//calculate 3vector magnitude:
double calculate3VecMag(double *vec) {
  double magnitude = sqrt(pow(vec[0],2.0) + pow(vec[1],2.0) + pow(vec[2],2.0));
  return magnitude;
}

//return the correctly scaled vector given its unit vector and magnitude:
void returnProperSizedVec(double *unitVec, double magnitude, double *vec) {
  vec[0] = unitVec[0]*magnitude;
  vec[1] = unitVec[1]*magnitude;
  vec[2] = unitVec[2]*magnitude;
}

//calculate vector direction between field point and the moving electron position:
void calculatePoyntingVec(double *fieldPoint, double *electronPos, double *vector) {
  vector[0] = fieldPoint[0] - electronPos[0];
  vector[1] = fieldPoint[1] - electronPos[1];
  vector[2] = fieldPoint[2] - electronPos[2];
}

//calculate linear 3 velocity direction from electron position: 
//(this must be normalised and then multiplied by the linear velocity
//calculated from kinetic energy to give the correct size of the vector).
/*void calculate3VelDirection(double *posVec, double pitchAngle, double *linearVelVec) {
  linearVelVec[0] = 1.0; //choose a value
  linearVelVec[1] = -linearVelVec[0]*(posVec[0]/posVec[1]); //y direction will scale according to x direction, using the product rule
  linearVelVec[2] = cos(pitchAngle); //z magnitude is fully dependent on pitch angle
}*/

//calculate acceleration vector from electron position:
/*void calculateAccelVec(double *posVec, double *accelVec) {
  accelVec[0] = -posVec[0];
  accelVec[1] = -posVec[1];
  accelVec[2] = -posVec[2];
}*/

//calculate the 4 separation:
void calculate4Separation(double R, double *observerPos, double *sourcePos, double *fourSeparation) {
  //double c = 299792458.0;
  fourSeparation[0] = R;
  fourSeparation[1] = observerPos[0] - sourcePos[0];
  fourSeparation[2] = observerPos[1] - sourcePos[1];
  fourSeparation[3] = observerPos[2] - sourcePos[2];
}

//calculate source's 4 velocity:
void calculate4Velocity(double lorentzFactor, double *threeVelocity, double *fourVelocity) {
  //double c = 299792458.0;
  fourVelocity[0] = lorentzFactor*c;
  fourVelocity[1] = lorentzFactor*threeVelocity[0];
  fourVelocity[2] = lorentzFactor*threeVelocity[1];
  fourVelocity[3] = lorentzFactor*threeVelocity[2];
}

//calculate the source's 4 acceleration:
void calculate4Accel(double lorentzFactor, double *threeVelocity, double *threeAccel, double *fourAccel) {
  //double c = 299792458.0;
  fourAccel[0] = pow(lorentzFactor,4.0)*dotProduct3vectors(threeAccel,threeVelocity)/c;
  fourAccel[1] = pow(lorentzFactor,2.0)*threeAccel[0] + 
                 (pow(lorentzFactor,4.0)*dotProduct3vectors(threeAccel,threeVelocity)/pow(c,2.0))*threeVelocity[0];
  fourAccel[2] = pow(lorentzFactor,2.0)*threeAccel[1] + 
                 (pow(lorentzFactor,4.0)*dotProduct3vectors(threeAccel,threeVelocity)/pow(c,2.0))*threeVelocity[1];;
  fourAccel[3] = pow(lorentzFactor,2.0)*threeAccel[2] + 
                 (pow(lorentzFactor,4.0)*dotProduct3vectors(threeAccel,threeVelocity)/pow(c,2.0))*threeVelocity[2];;
}

void Efield(double *fourSeparation, double *fourVelocity, double *fourAccel, 
            double R, double lorentzFactor, double *Efield) {
  //double c = 299792458.0;
  //double mu = 1.25663706e-6; //permeability of free space
  //double e = 1.602176634*pow(10.0,-19.0);
  double *fourVelSpatialComponents, *fourSeparationSpatialComponents;
  fourVelSpatialComponents = (double*)malloc(3*sizeof(double));
  fourSeparationSpatialComponents = (double*)malloc(3*sizeof(double));
  int i;
  for (i=0;i<3;i++) {
    fourVelSpatialComponents[i] = fourVelocity[i+1];
    fourSeparationSpatialComponents[i] = fourSeparation[i+1];
  }
  double Rs = lorentzFactor*(R - dotProduct3vectors(fourSeparationSpatialComponents,fourVelSpatialComponents)/c);
  double factor1 = (mu*e*pow(c,2.0))/(4*M_PI);
  double factor2 = (1.0 - dotProduct4vectors(fourSeparation,fourAccel)/pow(c,2.0))*(1.0/(c*pow(Rs,3.0)));
  double factor3 = 1.0/(c*pow(Rs,2.0));
    
  //Filling E field components:
  Efield[0] = factor1*(fourSeparation[1]*fourVelocity[0]-fourSeparation[0]*fourVelocity[1])*factor2 +
              (fourSeparation[1]*fourAccel[0]-fourSeparation[0]*fourAccel[1])*factor3;
  Efield[1] = factor1*(fourSeparation[2]*fourVelocity[0]-fourSeparation[0]*fourVelocity[2])*factor2 +
              (fourSeparation[2]*fourAccel[0]-fourSeparation[0]*fourAccel[2])*factor3;
  Efield[2] = factor1*(fourSeparation[3]*fourVelocity[0]-fourSeparation[0]*fourVelocity[3])*factor2 +
              (fourSeparation[3]*fourAccel[0]-fourSeparation[0]*fourAccel[3])*factor3;
  free(fourVelSpatialComponents); free(fourSeparationSpatialComponents);
}

void Bfield(double *fourSeparation, double *fourVelocity, double *fourAccel, 
            double R, double lorentzFactor, double *Bfield) {
  double *fourVelSpatialComponents, *fourSeparationSpatialComponents;
  fourVelSpatialComponents = (double*)malloc(3*sizeof(double));
  fourSeparationSpatialComponents = (double*)malloc(3*sizeof(double));
  int i;
  for (i=0;i<3;i++) {
    fourVelSpatialComponents[i] = fourVelocity[i+1];
    fourSeparationSpatialComponents[i] = fourSeparation[i+1];
  }
  double Rs = lorentzFactor*(R - dotProduct3vectors(fourSeparationSpatialComponents,fourVelSpatialComponents)/c);
  double factor1 = (mu*e*c)/(4*M_PI);
  double factor2 = (1.0 - dotProduct4vectors(fourSeparation,fourAccel)/pow(c,2.0))*(1.0/(c*pow(Rs,3.0)));
  double factor3 = 1.0/(c*pow(Rs,2.0));
    
  //Filling B field components:
  Bfield[0] = factor1*(fourSeparation[3]*fourVelocity[2]-fourSeparation[2]*fourVelocity[3])*factor2 +
              (fourSeparation[3]*fourAccel[2]-fourSeparation[2]*fourAccel[3])*factor3;
  Bfield[1] = factor1*(fourSeparation[1]*fourVelocity[3]-fourSeparation[3]*fourVelocity[1])*factor2 + 
              (fourSeparation[1]*fourAccel[3]-fourSeparation[3]*fourAccel[1])*factor3;
  Bfield[2] = factor1*(fourSeparation[2]*fourVelocity[1]-fourSeparation[1]*fourVelocity[2])*factor2 + 
              (fourSeparation[2]*fourAccel[1]-fourSeparation[1]*fourAccel[2])*factor3;

  free(fourVelSpatialComponents); free(fourSeparationSpatialComponents);
}



//magnitude of measured poynting vector for relativistic speeds:
/*double poyntingMagnitude(double *Efield, double *Bfield, double effArea) {
  double c = 299792458.0;
  double mu = 1.25663706e-6; //permeability of free space
  double EfieldMag = sqrt(pow(Efield[0],2.0) + pow(Efield[1],2.0) + pow(Efield[2],2.0));
  double BfieldMag = sqrt(pow(Bfield[0],2.0) + pow(Bfield[1],2.0) + pow(Bfield[2],2.0));
  double poyntingMag = (EfieldMag*BfieldMag)/(mu)*effArea;
  return poyntingMag;
}*/
double poyntingMagnitude(double *Efield, double effArea) {
  //double c = 299792458.0;
  //double mu = 1.25663706e-6; //permeability of free space
  double EfieldMag = sqrt(pow(Efield[0],2.0) + pow(Efield[1],2.0) + pow(Efield[2],2.0));
  double poyntingMag = (pow(EfieldMag,2.0))/(mu*c)*effArea;
  return poyntingMag;
}

int main() { 
  int i; 
  //double totalTime = 51*5.0e-5;
  double dt = 1.0e-9; //timestep
  int N = 500;
  //int N = 10000; //10,000 points
  //int N = 1200000; //Number of points to calculate
  double B = 1.0; //Idealised magnetic field
  double effArea = 1.0; //effective area of antenna, for idealised purposes set to 1
  //int N = totalTime/dt;
  double kineticEnergy = 18.575 * 1000.0 * e; //energy in Joules
  //double kineticEnergy =  100.0 * 1000.0 * e; //very relativistic electron
  //double kineticEnergy = 2500 * e; //beta 0.1 threshold is ~ 2578 electronVolts
  
  double *trajectoryPos, *trajectoryVel, *trajectoryAccel, *observerPos; //Input an electron trajectory, velocity and acceleration
                                       // and initialise observer position
  trajectoryPos = (double *)malloc(3*sizeof(double));
  trajectoryVel = (double *)malloc(3*sizeof(double));
  trajectoryAccel = (double *)malloc(3*sizeof(double));
  observerPos = (double *)malloc(3*sizeof(double));
  observerPos[0] = 0.0;
  observerPos[1] = 0.05;
  observerPos[2] = 0.0;

  //double pitchAngle = 87.0*(M_PI/180.0);
  double pitchAngle = 90.0*(M_PI/180.0);
   
  double *poyntingVec, *unitPoyntingVec, *unitVelVec, *unitAccelVec; 
  poyntingVec = (double *)malloc(3*sizeof(double));
  unitPoyntingVec = (double *)malloc(3*sizeof(double)); //unit vector pointing from source to observer, this is direction of our poynting vector.
  //Unit vectors to find direction of the 3velocity and 3 acceleration
  unitVelVec = (double *)malloc(3*sizeof(double));
  unitAccelVec = (double *)malloc(3*sizeof(double));


  double R; //distance of field point to electron --- will have to vary
  double *fourSeparation;
  double *threeVel, *threeAccel;
  double *fourVel, *fourAccel;
  fourSeparation = (double *)malloc(4*sizeof(double));
  threeVel = (double *)malloc(3*sizeof(double));
  threeAccel = (double *)malloc(3*sizeof(double));
  fourVel = (double *)malloc(4*sizeof(double));
  fourAccel = (double *)malloc(4*sizeof(double));
  
  double *EcalcField, *BcalcField;
  EcalcField = (double *)malloc(3*sizeof(double));
  BcalcField = (double *)malloc(3*sizeof(double));
  double *measuredPow, *phase;
  measuredPow = (double *)malloc(N*sizeof(double));
  phase = (double *)malloc(N*sizeof(double));
  double powerReceived;

  double **EfieldArray, **BfieldArray;
  EfieldArray = (double **)malloc(N*sizeof(double*));
  BfieldArray = (double **)malloc(N*sizeof(double*));
  for (i=0;i<N;i++) {
    EfieldArray[i] = (double *)malloc(3*sizeof(double));
    BfieldArray[i] = (double *)malloc(3*sizeof(double));
  }

  //determine whether to calculate beta from the relativistic kinetic energy, to choose
  //whether to use either the relativistic or non-relativistic equations:
  double beta = calculateBeta(kineticEnergy);
  double electronVel = beta*c;
  double lorentzFactor = 1.0/sqrt(1.0 - pow(beta,2.0));
  //If the electrons velocity is larger than 1/10th the speed of light, use the
  //relativistic calculation, otherwise, use the non-relativistic calculation:
  if (beta > 0.1) {
    printf("Using relativitic calculations...\n");
    printf("velocity calculated: %.4f\n", electronVel);
    double radius = (lorentzFactor*eMass*electronVel)/(e*B); //Lamor radius/gyroradius/cyclotron radius for relativistic motion
    double angularV = electronVel/radius;
    double accelMag = (e*electronVel*B)/(eMass*lorentzFactor);
    printf("lorentzFactor: %.4e\n", lorentzFactor);
    printf("radius: %.4e\n", radius);
    //read trajectory file:
    char fileName[40] = "idealTrajectoryDeltaT=";
    char fileType[5] = ".tsv";
    char cutoffVal[15];
    snprintf(cutoffVal, 15, "%.3e", dt);
    char *fileNameExtension;
    fileNameExtension = malloc(strlen(fileName)+1+20);
    strcpy(fileNameExtension, fileName);
    strcat(fileNameExtension, cutoffVal);
    strcat(fileNameExtension, fileType);
    printf("File read: %s", fileNameExtension);
    //char fileNameExtension[40] = "trajectoriesCyclotron.tsv";
    FILE *posData = fopen(fileNameExtension, "r");
    for (i=0;i<N;i++) {
      //Read data:
      //fscanf(posData, "%lf \t %lf \t %lf \n", &trajectoryPos[0], &trajectoryPos[1], &trajectoryPos[2]);
      //phase[i] = angularV*dt*i;
      fscanf(posData, "%lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \n",
                       &phase[i],
                       &trajectoryPos[0], &trajectoryPos[1], &trajectoryPos[2],
                       &trajectoryVel[0], &trajectoryVel[1], &trajectoryVel[2],
                       &trajectoryAccel[0], &trajectoryAccel[1], &trajectoryAccel[2]);
      //printf("position: %.3e \t %.3e \t %.3e\n", trajectoryPos[0], trajectoryPos[1], trajectoryPos[2]);

      //Find 3velocity from position:
      /*calculate3VelDirection(trajectoryPos, pitchAngle, threeVel); //finds 3velocity vector direction
      calculateUnitVec(threeVel, unitVelVec);  //normalises direction
      returnProperSizedVec(unitVelVec, electronVel, threeVel); //scales 3velocity vector correctly*/

      //Find 3acceleration from position:
      /*calculateAccelVec(trajectoryPos, threeAccel);  //finds 3acceleration vector direction
      calculateUnitVec(threeAccel, unitAccelVec);   //normalises direction
      returnProperSizedVec(unitAccelVec, accelMag, threeAccel); //scales 3acceleration vector correctly*/
      for (i=0;i<3;i++) {
        threeVel[i] = trajectoryVel[i];
      }

      for (i=0;i<3;i++) {
        threeAccel[i] = trajectoryAccel[i];
      }
      //printf("Dot product of a and v: %1.e\n", dotProduct3vectors(threeAccel,threeVel));
      
      //Calculate distance from source to observe point:
      R = calculateR(observerPos, trajectoryPos);
      
      //Calculate magnitude of position vector:
      //R = calculate3VecMag(trajectoryPos);

      //Calculate 4separation, X
      calculate4Separation(R, observerPos, trajectoryPos, fourSeparation); //<-fourSeparation = X
      //Calculate 4velocity, u
      calculate4Velocity(lorentzFactor, threeVel, fourVel); //<-fourVel = U
      //Calculate 4acceleration, a
      calculate4Accel(lorentzFactor, threeVel, threeAccel, fourAccel); //<-fourAccel = A

      //Calculate poynting vector direcion:
      calculatePoyntingVec(observerPos, trajectoryPos, poyntingVec);  //finds vector direction pointing to observer
      calculateUnitVec(poyntingVec, unitPoyntingVec);  //finds unit vector pointing to observer (for poynting vector)

      //Calculate Efield and Bfield from Faraday Tensor:
      Efield(fourSeparation, fourVel, fourAccel, R, lorentzFactor, EcalcField);
      Bfield(fourSeparation, fourVel, fourAccel, R, lorentzFactor, BcalcField);
      //Calculate poynting vector magnitude:
      //powerReceived = poyntingMagnitude(EcalcField, BcalcField, effArea);
      powerReceived = poyntingMagnitude(EcalcField, effArea);
      EfieldArray[i][0] = EcalcField[0];
      EfieldArray[i][1] = EcalcField[1];
      EfieldArray[i][2] = EcalcField[2];
      /*BfieldArray[i][0] = BcalcField[0];
      BfieldArray[i][1] = BcalcField[1];
      BfieldArray[i][2] = BcalcField[2];*/
      measuredPow[i] = powerReceived;
      if(i==N/10){
        printf("10 percent\n");
      }
      if(i==2*N/10){
        printf("20 percent\n");
      }
      if(i==3*N/10){
        printf("30 percent\n");
      }
      if(i==4*N/10){
        printf("40 percent\n");
      }
      if(i==N/2){
        printf("50 percent\n");
      }
      if(i==6*N/10){
        printf("60 percent\n");
      }
      if(i==7*N/10){
        printf("70 percent\n");
      }
      if(i==8*N/10){
        printf("80 percent\n");
      }
      if(i==9*N/10){
        printf("90 percent\n");
      }
    }
    printf("100 percent\n");
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
      fprintf(powData, "%.3e \t %.3e \n", phase[i], measuredPow[i]);
    }
    fclose(powData);
    free(fileNameExtension1);

    /*FILE *EMData = fopen("electricAndMagneticFields.tsv", "w");
    fprintf(EMData, "%s \t %s \t %s \t %s \t %s \t %s \n",
            "E_x", "E_y", "E_z", "B_x", "B_y", "B_z");
    for (i=0;i<N;i++) {
      fprintf(EMData, "%.4e \t %.4e \t %.4e \t %.4e \t %.4e \t %.4e \n ",
             EfieldArray[i][0], EfieldArray[i][1], EfieldArray[i][1], 
             BfieldArray[i][0], BfieldArray[i][1], BfieldArray[i][2]);
    }
    fclose(EMData);*/

    FILE *EData = fopen("electricAndMagneticFields.tsv", "w");
    fprintf(EData, "%s \t %s \t %s \n",
            "E_x", "E_y", "E_z");
    for (i=0;i<N;i++) {
      fprintf(EData, "%.4e \t %.4e \t %.4e \n ",
             EfieldArray[i][0], EfieldArray[i][1], EfieldArray[i][1]);
    }
    fclose(EData); 
  }


  for (i=0;i<N;i++) {
    free(EfieldArray[i]); free(BfieldArray[i]);
  }
  free(EfieldArray); free(BfieldArray); 
  free(threeVel); free(threeAccel);
  free(fourSeparation); free(fourVel); free(fourAccel);
  free(poyntingVec); free(unitPoyntingVec);
  free(unitVelVec); free(unitAccelVec);
  free(EcalcField); free(BcalcField);
  free(phase); free(measuredPow); free(trajectoryPos); free(trajectoryVel); free(trajectoryAccel); free(observerPos);
}