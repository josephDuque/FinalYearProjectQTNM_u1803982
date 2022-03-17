#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#define _USE_MATH_DEFINES
#include<math.h>

//Renamed "hertzianDipole.h" in the power measured code

//dot product of two 3-vectors:
double dotProduct3vectors(double *vec1, double *vec2){
  int i;
  double product = 0.0;
  for (i=0;i<3;i++){
    product = vec1[i]*vec2[i] + product;
  }
  return product;
}

//calculate 3vector magnitude:
double calculate3VecMag(double *vec) {
  double magnitude = sqrt(pow(vec[0],2.0) + pow(vec[1],2.0) + pow(vec[2],2.0));
  return magnitude;
}

//Calculate effective area factor from the separation vector, the wavelength and the antenna's orietation vector:
double calcEffArea(double *separation, double wavelength) {
  double *orientation;
  orientation = (double*)malloc(3*sizeof(double));
  orientation[0] = 0.0;
  orientation[1] = 0.0;
  orientation[2] = 1.0;
  double cosineVal = dotProduct3vectors(separation, orientation)/(calculate3VecMag(orientation)*calculate3VecMag(separation));  //Vectorial definition of the cosine
  double sinSquared = 1.0 - pow(cosineVal,2.0); //Trigonometric identity
  free(orientation);
  double effArea = (3/(8*M_PI))*pow(wavelength,2.0)*sinSquared;
  return effArea;
}