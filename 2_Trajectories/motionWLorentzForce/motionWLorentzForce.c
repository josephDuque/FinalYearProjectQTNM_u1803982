/* Calling DGELS using row-major order */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cblas.h>
#include <lapacke.h>
#include <string.h>
#include "magneticTrap.h"
//#include <sys/stat.h>
#define SIZE 3

//Function to calculate beta. The speed is then v=beta*c
double calculateBeta(double T) {
  double eMass = 9.10938356e-31;
  double c = 299792458.0;
  double factor = pow(T,2.0)/(pow(eMass,2.0)*pow(c,4.0)) + (2.0*T)/(eMass*pow(c,2.0)) + 1;
  double beta = sqrt(1.0 - 1.0/factor);
  return beta;
}

int main (int argc, const char * argv[]) {
    int i, j;
    //FILE *inputFile = fopen("simParms.dat","r");
    //for(int kk=0; kk<10; kk++){
    //}
    //
    double lambda,lambda0,Bx,By,Bz,q,dt,posX,posY,posZ;
    double kinE, mass,gamma,c,Ux,Uy,Uz,magV,magVComputed;
    double dirCosX, dirCosY, dirCosZ, rotationR, solenoidWidth;
    double coilRadius, current, mu, coilPosZMin1, coilPosZMax1, coilPosZMin2, coilPosZMax2;
    double Z_coil1_, Z_coil2_ ; 
    int    nTsteps, nCoils;

    //double vecRHS[SIZE]  = {0.0, 0.0, 0.0};
    //double Bfield[SIZE]  = {0.0, 0.0, 0.0 };
    //double Bcgrnd[SIZE]  = {0.0, 0.0, 1.0 };
    double *vecRHS, *Bfield, *background;
    vecRHS = (double*)malloc(SIZE*sizeof(double));
    Bfield = (double*)malloc(SIZE*sizeof(double));
    background = (double*)malloc(SIZE*sizeof(double));
    background[0] = 0.0; background[1] = 0.0; background[2] = 1.0;
    lapack_int info,nRowsA,nColsB,nColsARowsB,ldA,ldB,ipvt[SIZE*SIZE];

    dt= 5.0e-12;
    nTsteps = 2.5e6;
    q = 1.602176634e-19;
    c = 299792458.0;
    mass = 9.10938356e-31;
    kinE = 18575.0*q;
    rotationR = 4.6375e-4;
    solenoidWidth = 0.2;
    coilRadius = 0.05; 
    current = 20.0;
    nCoils = 11;
    mu = 1.25663706e-6;
    Z_coil1_ = -0.5;
    Z_coil2_ =  0.5;
    coilPosZMin1 = Z_coil1_ - solenoidWidth/2.0;
    coilPosZMax1 = Z_coil1_ + solenoidWidth/2.0;
    coilPosZMin2 = Z_coil2_ - solenoidWidth/2.0;
    coilPosZMax2 = Z_coil2_ + solenoidWidth/2.0;

    posX=rotationR;posY=0.0;posZ=0.0;
    bathTubFieldSolenoids(coilRadius, current, coilPosZMin1, coilPosZMax1, coilPosZMin2, coilPosZMax2, nCoils,
                           mu, posX, posY, posZ, background, Bfield);
    Bx = Bfield[0]; By = Bfield[1]; Bz = Bfield[2]; 
    //-----------
    magV = calculateBeta(kinE)*c;
    double pitchAngle = 85.0 * M_PI/180.0;
    dirCosX = 0.0; dirCosZ=cos(pitchAngle); dirCosY=1.0-dirCosZ;
    Ux = magV*dirCosX; Uy=magV*dirCosY; Uz=magV*dirCosZ;
    //double ufuture[SIZE] = {0.0, 0.0, 0.0};
	//double upast[SIZE]   = {Ux, Uy, Uz};
    double *ufuture, *upast;
    ufuture = (double*)malloc(SIZE*sizeof(double));
    upast = (double*)malloc(SIZE*sizeof(double));
    upast[0] = Ux; upast[1] = Uy; upast[2] = Uz;
    printf("# Starting simulation with U(x), U(y), U(z) = %8.5e , %8.5e , %8.5e \n", Ux, Uy, Uz);
    gamma = 1.0/(sqrt(1.0-(pow(magV,2.0))/(pow(c,2.0))));
    lambda0= (0.5*q*dt)/(2.0*mass*gamma); // Careful: this definition is only for the first time step
    lambda = (1.0*q*dt)/(2.0*mass*gamma); // Careful: this definition is only for the first time step

    //writing a and b as vectors, but they are in concept the matrices of the system of linear equations
    //double a[SIZE*SIZE]  = { 1.0, lambda0*Bz,-lambda0*By , -lambda0*Bz,1.0, lambda0*Bx ,  lambda0*By,-lambda0*Bx,1.0 };
    //double b[SIZE*SIZE]  = { 1.0,-lambda0*Bz, lambda0*By ,  lambda0*Bz,1.0,-lambda0*Bx , -lambda0*By, lambda0*Bx,1.0 };
    double *a, *b;
    a = (double*)malloc(SIZE*SIZE*sizeof(double));
    b = (double*)malloc(SIZE*SIZE*sizeof(double));

    a[0] = 1.0;         a[1] = lambda0*Bz;  a[2] = -lambda0*By;
    a[3] = -lambda0*Bz; a[4] = 1.0;         a[5] = lambda0*Bx;
    a[6] = lambda0*By;  a[7] = -lambda0*Bx; a[8] = 1.0; 

    b[0] = 1.0;         b[1] = -lambda0*Bz; b[2] = lambda0*By;
    b[3] = lambda0*Bz;  b[4] = 1.0;         b[5] = -lambda0*Bx;
    b[6] = -lambda0*By; b[7] = lambda0*Bx;  b[8] = 1.0;


    nRowsA = SIZE;
    nColsB = 1;
    nColsARowsB = SIZE;
    ldA    = nRowsA;
    ldB    = nColsARowsB;

    // Printing LHS matrix
    /*printf("#LHS Matrix (Traspose -> Column major order -----------\t | \n");
    for(i=0;i<SIZE;i++)
    {
        printf("#| \t");
        for(j=0;j<SIZE;j++)
        {
           printf("%.5e \t",a[i*SIZE + j]);
        }
        printf("\t | \n");
    }
    printf("#------------------------------------------------------\n");
    printf("#RHS Matrix (Traspose -> Column major order -----------\t | \n");
    for(i=0;i<SIZE;i++)
    {
        printf("#| \t");
        for(j=0;j<SIZE;j++)
        {
           printf("%.5e \t",b[i*SIZE + j]);
        }
        printf("\t | \n");
    }
    printf("#----------------------------------------------------------\n");
    printf("# Time       | Ux \t        | Uy \t        | Uz \t        | magV \t        | PosX \t        | PosY \t        | PosZ \t        |\n");
    */
    double timeVal = 0.0;
    for (int tt = 0; tt < nTsteps; tt++) {
        timeVal = (double)tt * dt;
        //using header variables from LAPACK and openblas
        cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
            nRowsA,nColsB,nColsARowsB,1.0,b,ldA,upast,ldB,0.0,vecRHS,nRowsA);
        for (i=0; i<SIZE; i++){ufuture[i] = vecRHS[i];}; // This is because dgesv rewrites ufuture
        info = LAPACKE_dgesv(LAPACK_COL_MAJOR, SIZE, 1, a, ldA, ipvt, ufuture, ldA);
        for (i=0; i<SIZE; i++){upast[i] = ufuture[i];};
        Ux = ufuture[0]; 
        Uy = ufuture[1];
        Uz = ufuture[2];
        magVComputed = sqrt(pow(Ux,2.0)+pow(Uy,2.0)+pow(Uz,2.0));
        magV = magVComputed;
        gamma = 1.0/(sqrt(1.0-(pow(magV,2.0))/(pow(c,2.0))));
        lambda = (q*dt)/(2.0*mass*gamma);

        bathTubFieldSolenoids(coilRadius, current, coilPosZMin1, coilPosZMax1, 
                coilPosZMin2, coilPosZMax2, nCoils, mu, posX, posY, posZ, background, Bfield);                
        Bx = Bfield[0]; By = Bfield[1]; Bz = Bfield[2]; 

        a[0] = 1.0;        a[1] = lambda*Bz;    a[2] = -lambda*By; 
        a[3] = -lambda*Bz; a[4] = 1.0;          a[5] = lambda*Bx;
        a[6] = lambda*By;  a[7] = -lambda*Bx;   a[8] = 1.0;

        b[0] = 1.0;        b[1] = -lambda*Bz;   b[2] = lambda*By; 
        b[3] = lambda*Bz;  b[4] = 1.0;          b[5] = -lambda*Bx;
        b[6] = -lambda*By; b[7] = lambda*Bx;    b[8] = 1.0;
        // -----------
        posX = posX + dt*Ux;
        posY = posY + dt*Uy;
        posZ = posZ + dt*Uz;
        if ( tt%1 == 0 ) {    //This alters the number of points written to the text file
            printf(" %7.4e \t %8.5e \t %8.5e \t %8.5e \t %8.5e \t %8.5e \t %8.5e \n ",
                timeVal, posX, posY, posZ, ufuture[0], ufuture[1], ufuture[2]);
        }
    }

    //free(a); free(b);
    free(ufuture); free(upast);
    free(vecRHS); free(Bfield); free(background);
    info = 0;
    return(info);
}
