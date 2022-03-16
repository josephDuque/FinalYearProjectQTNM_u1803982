#include <stdio.h>
#include <stdlib.h>
#define _use_math_defines
#include <math.h>

double centralField(double radius, double current, double mu) {
  double centralField = (current*mu)/(2.0*radius);
  return centralField;
}

double onAxisField(double radius, double current, double Z, double mu, double z_p) {
  double onAxisField = (mu*current*pow(radius,2.0))/(2.0*(pow(radius,2.0) + pow((z_p - Z),2.0),1.5));
  return onAxisField;
}

double ellipIntegrandK(double beta, double t) {
  double integrand = pow((1.0 - pow((beta*sin(t)),2.0)),-0.5);
  return integrand;
}

double ellipIntegrandE(double beta, double t) {
  double integrand = pow((1.0 - pow((beta*sin(t)),2.0)),0.5);
  return integrand;
}

double variableChange(double x) {
  double t = (M_PI/4)*(x + 1.0);
  return t;
}

void cirularCoilField(double radius, double current, double Z, double mu, 
                      double x_p, double y_p, double z_p, double *Bfield) {
  double R = sqrt(pow(x_p,2.0)+pow(y_p,2.0)); //radius of point
  //If point is on axis:
  if (R/radius <1e-10) {
    Bfield[0] = 0.0; //x coordinate of B field
    Bfield[1] = 0.0; //y coordinate of B field
    Bfield[2] = onAxisField(radius, current, Z, mu, z_p); //z coordinate of B field
  }
  //If not on axis:
  double z_rel = z_p - Z; //z relative to position of coil
  double b_central = centralField(radius, current, mu);
  //normalised radius of point and normalised z position.
  double r_hat = R/radius;
  double z_hat = z_rel/radius;
  
  double alpha = pow((1.0 + r_hat),2.0) + pow(z_hat,2.0);
  double piRootAlpha = M_PI*sqrt(alpha);
  double beta = (4.0*r_hat)/alpha;
  double gamma = alpha - 4.0*r_hat;

  double ellipIntK, ellipIntE;
  //finding elliptic integrals using gaussian quadrature:
  double t1 = variableChange(1.0/sqrt(3.0));
  double t2 = variableChange(-1.0/sqrt(3.0));
  //two point quadrature:
  ellipIntK = (M_PI/4.0)*(ellipIntegrandK(beta,t1) 
               + ellipIntegrandK(beta,t2));  
  ellipIntE = (M_PI/4.0)*(ellipIntegrandE(beta,t1)
               + ellipIntegrandE(beta,t2));  
  
  double b_r = b_central * (ellipIntE * ((1.0 + pow(r_hat,2.0) + pow(z_hat,2.0))/gamma) - ellipIntK)/piRootAlpha * (z_rel/R);
  double b_z = b_central * (ellipIntE * ((1.0 - pow(r_hat,2.0) - pow(r_hat,2.0))/gamma) + ellipIntK)/piRootAlpha;
  
  Bfield[0] = (b_r*x_p)/R;
  Bfield[1] = (b_r*y_p)/R;
  Bfield[2] = b_z;
}

void solenoidCoilField(double radius, double current, double Zmin, double Zmax, int Ncoils,
                       double mu, double x_p, double y_p, double z_p, double *Bfield) {
  double *tempCoil, *tempSolenoid;
  int i,j;
  tempCoil = (double *)malloc(3*sizeof(double));
  tempSolenoid = (double *)malloc(3*sizeof(double));
  for (i=0;i<3;i++) {
    tempSolenoid[i] = 0.0;
  }
  double step = (Zmax-Zmin)/(Ncoils-1);
  double ZcoilPos;
  for (i=0;i<Ncoils;i++) {
    ZcoilPos = Zmin+(i*step);
    cirularCoilField(radius, current, ZcoilPos, mu, x_p, y_p, z_p, tempCoil);
    for (j=0;j<3;j++) {
      tempSolenoid[j] = tempSolenoid[j] + tempCoil[j];
    }
  }
  for (i=0;i<3;i++) {
    Bfield[i] = tempSolenoid[i];
  }
  free(tempCoil); free(tempSolenoid);
}

//bathtub trap of one coil on either side:
void bathTubField1Coil(double radius, double current, double z_coil1, double z_coil2,
                  double mu, double x_p, double y_p, double z_p,
		  double *background, double *bathField) {
  double *tempCoil1, *tempCoil2;
  int i;
  tempCoil1 = (double *)malloc(3*sizeof(double));
  tempCoil2 = (double *)malloc(3*sizeof(double));
  cirularCoilField(radius, current, z_coil1, mu, x_p, y_p, z_p, tempCoil1);
  cirularCoilField(radius, current, z_coil2, mu, x_p, y_p, z_p, tempCoil2);
  
  for (i=0;i<3;i++) {
    bathField[i] = tempCoil1[i] + tempCoil2[i] + background[i];
  }
  free(tempCoil1); free(tempCoil2);
}

//bathtub trap for solenoid coils on either side:
void bathTubFieldSolenoids(double radius, double current, double z_coil1min, double z_coil1max,
                  double z_coil2min, double z_coil2max, int Ncoils,
		  double mu, double x_p, double y_p, double z_p,
		  double *background, double *bathField) {
  double *tempCoil1, *tempCoil2;
  int i;
  tempCoil1 = (double *)malloc(3*sizeof(double));
  tempCoil2 = (double *)malloc(3*sizeof(double));
  solenoidCoilField(radius, current, z_coil1min, z_coil1max, Ncoils, mu, x_p, y_p, z_p, tempCoil1);
  solenoidCoilField(radius, current, z_coil2min, z_coil2max, Ncoils, mu, x_p, y_p, z_p, tempCoil2);

  for (i=0;i<3;i++) {
    bathField[i] = tempCoil1[i] + tempCoil2[i] + background[i];
  }
  free(tempCoil1); free(tempCoil2);
}

int main() {
  int i;
  double *background;
  background = (double *)malloc(3*sizeof(double));
  //1 Tesla background field
  for (i=0;i<3;i++) {
    background[i] = 1.0;
  }
  double radius = 0.005;
  double current = 40.0;
  double Z_coil1 = -1.0;
  double Z_coil2 = 1.0;
  double mu = 1.25663706e-6; //permeability of free space
  
  //point in space to probe:
  double xPoint;
  double yPoint;  
  double zPoint;
  //calculating single coil bathtub field on the z axis:--------------------------------
  xPoint = 0;
  yPoint = 0;
  int N = 1001;
  double *zArray, *bathField, *BzArray;
  bathField = (double *)malloc(3*sizeof(double));
  zArray = (double *)malloc(N*sizeof(double));
  BzArray = (double *)malloc(N*sizeof(double));
  double trapLength = 4.0;
  double zStep = trapLength/(N-1);
  for (i=0;i<N;i++) {
    zArray[i] = -2.0+(i*zStep);
    bathTubField1Coil(radius, current, Z_coil1, Z_coil2, mu, xPoint, yPoint, zArray[i], background, bathField);
    BzArray[i] = bathField[2];
  }
  //printing:
  FILE *bathtubZ = fopen("bathtubZaxis_1coil.tsv", "w");
  for (i=0;i<N;i++) {
    fprintf(bathtubZ, "%f \t %f \n", zArray[i], BzArray[i]);
  }
  fclose(bathtubZ);
  free(zArray); free(BzArray); free(bathField);
  //free(zArray); 
  //------------------------------------------------------------------------------------
  //calculating solenoid coil bathtub field on the z axis:------------------------------
  //xPoint = 0;
  //yPoint = 0;
  //int N = 1001;
  double *zArray2, *BzArray2, *bathField2;
  bathField2 = (double *)malloc(3*sizeof(double));
  zArray2 = (double *)malloc(N*sizeof(double));
  BzArray2 = (double *)malloc(N*sizeof(double));
  int Ncoils = 11;
  double solenoidWidth = 0.2;
  double Z_coil1min = Z_coil1 - solenoidWidth/2.0;
  double Z_coil1max = Z_coil1 + solenoidWidth/2.0;
  double Z_coil2min = Z_coil2 - solenoidWidth/2.0;
  double Z_coil2max = Z_coil2 + solenoidWidth/2.0;
  //double trapLength = 4.0;
  //double zStep = trapLength/(N-1);
  for (i=0;i<N;i++) {
    zArray2[i] = -2.0+(i*zStep);
    bathTubFieldSolenoids(radius, current, Z_coil1min, Z_coil1max, Z_coil2min, Z_coil2max, Ncoils,
                           mu, xPoint, yPoint, zArray2[i], background, bathField2);
    BzArray2[i] = bathField2[2];
  }
  //printing:
  FILE *bathtubZsol = fopen("bathtubZaxis_solenoidCoil.tsv", "w");
  for (i=0;i<N;i++) {
    fprintf(bathtubZsol, "%f \t %f \n", zArray2[i], BzArray2[i]);
  }
  fclose(bathtubZsol);
  free(zArray2); free(BzArray2); free(bathField2);
  //------------------------------------------------------------------------------------
  free(background); 
}