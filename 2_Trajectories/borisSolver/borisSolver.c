#include <stdio.h>
#include <stdlib.h>
#define _use_math_defines
#include <math.h>
#include "magneticTrap.h"

//Calculate fraction of speed of light:
double calculateBeta(double T) {
  double eMass = 9.10938356*pow(10.0,-31.0);
  double c = 299792458.0;

  double factor = pow(T,2.0)/(pow(eMass,2.0)*pow(c,4.0)) + (2.0*T)/(eMass*pow(c,2.0)) + 1;
  double beta = sqrt(1.0 - 1.0/factor);
  return beta;
}

//Finds the cyclotron frequency in component form, from position:
void calcOmega(double T, double *B, double *omega) {
  double e = 1.602176634*pow(10.0,-19.0);
  double eMass = 9.10938356*pow(10.0,-31.0);
  double c = 299792458.0;
  int i;
  for (i=0;i<3;i++) {
    omega[i] = ((e * B[i])/(eMass*pow(10.0,3.0) + T/pow(c,2)))/(2.0*M_PI);
  }
}

//calculate 3vector magnitude:
double calculate3VecMag(double *vec) {
  double magnitude = sqrt(pow(vec[0],2.0) + pow(vec[1],2.0) + pow(vec[2],2.0));
  return magnitude;
}

// cross product of two 3-vectors:
void crossProduct3vector(double *vec1, double *vec2, double *cross_P) {
    cross_P[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
    cross_P[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
    cross_P[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];
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

//acceleration from radiation reaction force: 
//Uses Ford O'Connel equation in 3D
void radiationAccel(double tau, double *omega,
		    double *vel, double *accel) {

   double denom = 1.0/(1.0 + pow(tau,2.0) * dotProduct3vectors(omega, omega));
   accel[0] = 0.0;
   accel[1] = 0.0;
   accel[2] = 0.0;

   //Larmor terms for acceleration:
   accel[0] = accel[0] - tau * (pow(omega[2],2.0) + pow(omega[1],2.0)) * vel[0];
   accel[0] = accel[0] + tau * omega[0] * (omega[2]*vel[2] + omega[1]*vel[1]);
   accel[1] = accel[1] - tau * (pow(omega[2],2.0) + pow(omega[0],2.0)) * vel[1];
   accel[1] = accel[1] + tau * omega[1] * (omega[2]*vel[2] + omega[0]*vel[0]);
   accel[2] = accel[2] - tau * (pow(omega[0],2.0) + pow(omega[1],2.0)) * vel[2];
   accel[2] = accel[2] + tau * omega[2] * (omega[0]*vel[0] + omega[1]*vel[1]);

   accel[0] = accel[0]*denom;
   accel[1] = accel[1]*denom;
   accel[2] = accel[2]*denom;
}
//-------------------------------------------------------------------------------------//
//-------------------------------------------------------------------------------------//
//-------------------------------------------------------------------------------------//
//-------------------------------------------------------------------------------------//
//-------------------------------------------------------------------------------------//
//Use Boris solver to advance equations of motion:
void advanceStep(double radius, double current, double z_coil1min, double z_coil1max,
                 double z_coil2min, double z_coil2max, int Ncoils,
		 double mu, double T,
		 double *background,
		 double tau, 
		 double time_step, double *initialPos, double *initialVel,
		 double *nextPos, double *nextVel) {
  double c = 299792458.0;
  double e = 1.602176634*pow(10.0,-19.0);
  double eMass = 9.10938356*pow(10.0,-31.0);
  int i;
  double gamma_n = 1.0/sqrt(1.0 - dotProduct3vectors(initialVel, initialVel)/pow(c,2.0));
  double *u_n, *x_n, *v_n, *omega;
  u_n = (double*)malloc(3*sizeof(double));    
  x_n = (double*)malloc(3*sizeof(double));    
  v_n = (double*)malloc(3*sizeof(double));    
  omega = (double*)malloc(3*sizeof(double));    
  for (i=0;i<3;i++) {
    u_n[i] = initialVel[i]*gamma_n;
    x_n[i] = initialPos[i];
    v_n[i] = initialVel[i];
  }

  //Half position step:
  double *x_nplushalf, *u_minus;
  x_nplushalf = (double*)malloc(3*sizeof(double)); 
  u_minus = (double*)malloc(3*sizeof(double));    
  for (i=0;i<3;i++) {
    x_nplushalf[i] = x_n[i] + v_n[i]*time_step/2.0;  //Half position
  }

  //Calculate B field at half-position
  double *B_nplushalf;
  B_nplushalf = (double*)malloc(3*sizeof(double));
  bathTubFieldSolenoids(radius, current, z_coil1min, z_coil1max,
                    		      z_coil2min, z_coil2max, Ncoils,
		                      mu, x_nplushalf[0], x_nplushalf[1], x_nplushalf[2],
		                      background, B_nplushalf);  
  //Get omega from half position:
  calcOmega(T, B_nplushalf, omega);
  //First half of the radiation reaction acceleration:
  double *accel;
  accel = (double*)malloc(3*sizeof(double)); 
  radiationAccel(tau, omega, 
		 u_n, accel);   //Use updated omega to calculate acceleration
  for (i=0;i<3;i++) {
    u_minus[i] = u_n[i] + (time_step/2.0)*accel[i];  //Calculate half position velocity
  }
  double gamma_minus = sqrt(1.0 + dotProduct3vectors(u_n,u_n)/pow(c,2.0));

  //Rotation step:
  double BMag = calculate3VecMag(B_nplushalf);
  double theta = (e*time_step)/(eMass*gamma_minus) * calculate3VecMag(B_nplushalf);
  double *t, *u_prime, *u_plus;
  t = (double*)malloc(3*sizeof(double));
  u_prime = (double*)malloc(3*sizeof(double));
  u_plus = (double*)malloc(3*sizeof(double));
  
  for (i=0;i<3;i++) { 
   t[i] = tan(theta/2.0) * B_nplushalf[i]/BMag;
  }
  double *u_minusCrossT, *u_primeCrossT;
  u_minusCrossT = (double*)malloc(3*sizeof(double));
  u_primeCrossT = (double*)malloc(3*sizeof(double));
  crossProduct3vector(u_minus, t, u_minusCrossT);
  for (i=0;i<3;i++) { 
    u_prime[i] = u_minus[i] + u_minusCrossT[i];
  }
  crossProduct3vector(u_prime, t, u_primeCrossT);
  for (i=0;i<3;i++) { 
    u_plus[i] = u_minus[i] + 2.0/(1.0 + dotProduct3vectors(t, t)) * u_primeCrossT[i];
  }
  double *u_nplus1, *v_nplus1, *x_nplus1; 
  u_nplus1 = (double*)malloc(3*sizeof(double)); 
  v_nplus1 = (double*)malloc(3*sizeof(double));
  x_nplus1 = (double*)malloc(3*sizeof(double));

  //Second half of the radiation reaction acceleration:
  radiationAccel(tau, omega, 
   	  	 u_plus, accel); //Position is not updated, so omega is not updated
  for (i=0;i<3;i++) { 
    u_nplus1[i] = u_plus[i] + (time_step/2.0) * accel[i];
  }
  double gamma_nplus1 = sqrt(1.0 + pow((calculate3VecMag(u_nplus1)/c),2.0));
  //Update position
  for (i=0;i<3;i++) { 
    nextVel[i] = u_nplus1[i] / gamma_nplus1;
    nextPos[i] = x_nplushalf[i] + nextVel[i] * time_step/2.0;
  }
  
  free(u_nplus1); free(v_nplus1); free(x_nplus1);
  free(u_minusCrossT); free(u_primeCrossT);
  free(t); free(u_prime); free(u_plus);
  free(x_nplushalf); free(u_minus); free(B_nplushalf);
  free(u_n); free(v_n); free(x_n); free(accel); free(omega);
}
//-------------------------------------------------------------------------------------//
//-------------------------------------------------------------------------------------//
//-------------------------------------------------------------------------------------//
//-------------------------------------------------------------------------------------//
//-------------------------------------------------------------------------------------//
//Solve equations using Boris method, for n rotations:
void solve(double radius, double current, 
	   double z_coil1min, double z_coil1max, double z_coil2min, double z_coil2max, int Ncoils,
	   double mu, double T,
	   double *background,
   	   double tau, double *omega, 
	   double time_step, int nSteps, double *initialPos, double *initialVel,
	   double *times, double **pos, double **vel,int n_rotations)  {
  double c = 299792458.0;

  /*double *u_n, *x_n, *v_n;
  u_n = v0*1.0 /sqrt(1.0 - dotProduct3vectors(v0,v0)/pow(c,2.0));
  x_n = x0;
  v_n = v0;*/
  
  pos[0][0] = initialPos[0];
  pos[1][0] = initialPos[1];
  pos[2][0] = initialPos[2];
  vel[0][0] = initialVel[0];
  vel[1][0] = initialVel[1];
  vel[2][0] = initialVel[2];

  //Loop through steps, using the advance_step function at each point
  int i;
  double *pos_n, *vel_n; //nth position and velocity at each timestep
  pos_n = (double*)malloc(3*sizeof(double));
  vel_n = (double*)malloc(3*sizeof(double));
  
  double *pos_nminus1, *vel_nminus1;
  pos_nminus1 = (double*)malloc(3*sizeof(double));
  vel_nminus1 = (double*)malloc(3*sizeof(double));
  
  for (i=1;i<nSteps;i++) {
    pos_nminus1[0] = pos[0][i-1]; /**/ pos_nminus1[1] = pos[1][i-1]; /**/  pos_nminus1[2] = pos[2][i-1];
    vel_nminus1[0] = vel[0][i-1]; /**/ vel_nminus1[1] = vel[1][i-1]; /**/  vel_nminus1[2] = vel[2][i-1];
    advanceStep(radius, current, 
    		z_coil1min, z_coil1max, z_coil2min, z_coil2max, Ncoils,
		mu, T,
		background,
		tau,
		time_step, pos_nminus1, vel_nminus1,
		pos_n, vel_n); //Outputs the next position and velocity
    times[i] = i*time_step;
    pos[0][i] = pos_n[0];
    pos[1][i] = pos_n[1];
    pos[2][i] = pos_n[2];
    vel[0][i] = vel_n[0];
    vel[1][i] = vel_n[1];
    vel[2][i] = vel_n[2];
  }
  free(pos_nminus1); free(vel_nminus1);
  free(pos_n); free(vel_n);
}
//-------------------------------------------------------------------------------------//
//-------------------------------------------------------------------------------------//
//-------------------------------------------------------------------------------------//
//-------------------------------------------------------------------------------------//
//-------------------------------------------------------------------------------------//
//Calculate acceleration at a position due to
//Lorentz force and Larmor radiation reaction terms:
void calcAccel(double tau, double *omega, 
	       double *velocity, double *accel){
  //Lorentz force perpendicular to velocity and angular velocity:
  crossProduct3vector(velocity,omega,accel);
  int i;
  double *temp;
  temp = (double*)malloc(3*sizeof(double));
  //Larmor radiation reaction terms:
  radiationAccel(tau, omega,
   		 velocity, temp);
  for (i=0;i<3;i++) {
    accel[i] = accel[i] + temp[i];
  }
  free(temp);
}



int main() {
  int i;
  
  double mu = 1.25663706e-6; //permeability of free space
  double c = 299792458.0;
  double e = 1.602176634*pow(10.0,-19.0);
  
  double kineticEnergy = 18.575 * 1000 * e;

  //1 Tesla background field
  double *background;
  background = (double *)malloc(3*sizeof(double));
  background[0] = 0.0;  background[1] = 0.0;  background[2] = 1.0;
  
  double tau = 0.0; //Larmor power parameter

  double radius = 0.005; //radius of solenoids
  double current = 0.0; //If current=0, we just have a background field with no solenoids
  double Z_coil1 = -1.0;
  double Z_coil2 = 1.0;

  int Ncoils = 101;
  double solenoidWidth = 0.2;
  double z_coil1min = Z_coil1 - solenoidWidth/2.0;
  double z_coil1max = Z_coil1 + solenoidWidth/2.0;
  double z_coil2min = Z_coil2 - solenoidWidth/2.0;
  double z_coil2max = Z_coil2 + solenoidWidth/2.0;


  //initial position and inital velocity
  double *x0, *v0;
  x0 = (double*)malloc(3*sizeof(double));
  v0 = (double*)malloc(3*sizeof(double));

  x0[0] = 1.0; x0[1] = 0.0; x0[2] = 0.0; 
  v0[0] = 0.0; v0[1] = 1.0; v0[2] = 0.0; //This will be multiplied by magnitude of velocity from kinetic energy
  
  //Calculate Bfield at the starting positions to calculate the
  //cyclotron frequency at the start, omega
  double *Bfield;
  Bfield = (double*)malloc(3*sizeof(double));
  bathTubFieldSolenoids(radius, current, z_coil1min, z_coil1max,
                    	z_coil2min, z_coil2max, Ncoils,
		        mu, x0[0], x0[1], x0[2],
		        background, Bfield);
  printf("B: %.4e \t %.4e \t  %.4e\n", Bfield[0], Bfield[1], Bfield[2]);
  //Get magnitude of omega:
  double *omega;
  omega = (double*)malloc(3*sizeof(double)); 
  calcOmega(kineticEnergy, Bfield, omega);
  //printf("omega: %.4e \t %.4e \t  %.4e\n", omega[0], omega[1], omega[2]);
  double omega0 = calculate3VecMag(omega);
  //Correct omega magnitude to consider relativity:
  double lorentzFactor = 1.0/sqrt(1.0 - dotProduct3vectors(v0,v0)/pow(c,2.0));
  omega0 = omega0/lorentzFactor;

  //=========================================================/
  //Calculate time step:
  double factor=1.0e-3;  //<---Change this for dt
  double max_step = 2.0*M_PI*factor/omega0; //This is the time step
  //printf("omega0: %.4e\n", omega0);
  
  //Number of gyro-orbits { as function of B(position) }
  int n_rotations = 20; //<----Change this for N
  //=========================================================/

  //Final time:
  double t_end = n_rotations*2.0*M_PI/omega0;
  //Calculate number of steps:
  int nSteps = (int)ceil(t_end/max_step) + 1;
  printf("Number of steps: %d\n", nSteps);
  
  double time_step = max_step;
  printf("Time step: %.4e\n", time_step);

  //Allocating printing arrays:
  double *times;
  times = (double*)malloc(nSteps*sizeof(double));
  double **positions, **velocities;
  positions = (double**)malloc(nSteps*sizeof(double*));
  velocities = (double**)malloc(nSteps*sizeof(double*));
  for (i=0;i<nSteps;i++) {
    positions[i] = (double*)malloc(3*sizeof(double));
    velocities[i] = (double*)malloc(3*sizeof(double));
  }
  double *accel;
  accel = (double*)malloc(3*sizeof(double));
  
  //Performing Boris solver:
  solve(radius, current, 
	z_coil1min, z_coil1max, z_coil2min, z_coil2max, Ncoils,
	mu, kineticEnergy,
	background,
   	tau, omega, 
	time_step, nSteps, x0, v0,
	times, positions, velocities, n_rotations);

  double *currentPos, *currentVel;
  currentPos = (double*)malloc(3*sizeof(double));
  currentVel = (double*)malloc(3*sizeof(double));
  //printing time, positions, velocities, and accelerations:
  FILE *borisTrajectories = fopen("trajectoriesBoris.tsv", "w");
  for (i=0;i<nSteps;i++) {
    currentPos[0] = positions[i][0]; currentPos[1] = positions[i][1]; currentPos[2] = positions[i][2];
    currentVel[0] = velocities[i][0]; currentVel[0] = velocities[i][1]; currentVel[2] = velocities[i][2];	  
    bathTubFieldSolenoids(radius, current, z_coil1min, z_coil1max,
                    	  z_coil2min, z_coil2max, Ncoils,
		          mu, currentPos[0], currentPos[1], currentPos[2],
		          background, Bfield); //Calculate B field to update omega
    calcOmega(kineticEnergy, Bfield, omega); //Calculate omega at new position
    calcAccel(tau, omega, currentVel, accel); //Calculate acceleration at position
    fprintf(borisTrajectories, 
    	    "%.8e \t %.8e \t %.8e \t %.8e \t %.8e \t %.8e \t %.8e \t %.8e \t %.8e \t %.8e \n",
	     times[i], positions[i][0], positions[i][1], positions[i][2], 
	               velocities[i][0], velocities[i][1], velocities[i][2], 
		       accel[0], accel[1], accel[2]);
  }
  fclose(borisTrajectories);

  //Free arrays
  for (i=0;i<nSteps;i++){
    free(positions[i]);
    free(velocities[i]);
  }
  free(positions); free(velocities);
  free(currentPos); free(currentVel);
  free(accel);
  free(x0); free(v0);
  free(omega); free(background); free(Bfield);
}