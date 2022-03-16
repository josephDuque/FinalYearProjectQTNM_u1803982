# FinalYearProjectQTNM_u1803982
Repository with the code used to produce numerical result of the final year
project. The code included simulates the motion of an electron with a kinetic
energy of 18.575 keV inside a magnetic bathtub trap. The trajectories produced
are then inputted into code that calculates the power emitted at each position
in the trajectory, and saves these values to a text file against time and phase.
With the signals produced, data analysis code such as fourier transform code in
Matlab was employed. Thecode is split into modular pieces in order to make
aquisition of results and demonstration of the workflow easier.

The directory structure is as follows:

1_MagneticTrap: -Code written to produce the magnetic field in a magnetic trap
for a given position

2_Trajectories: Contains multiple folders:

  borisSolver: -code to solve the equation of motion as dicated by the Lorentz
  and Abaraham-Lorentz force. This code is only half finished and requires more
  polish to produce results.
  
  circularMotionEnergyLoss: -simple code to simulate a circular motion
  trajectory where the electron losses constant energy over time (constant power
  radiated assumption.
  
  idealCircularMotion: -simple code to produce circular motion trajectories.

  idealCyclotronMotion: -simple code superimposing circular motion with non-zero
  velocity in the z-axis. This code was used to produce a figure in the
  introduction of the report.

  motionWLorentzForce: -this is the principal code that returned the
  trajectories of an electron in a magnetic field. The mathematematical
  formulation of the algorithm implies solving a system of linear equations for
  each timestep. For this reason the code employed the use of the LAPACK library
  in order to solve the system with the use of linear algebra methods. 
  This code contains a Makefile in order to compile and run the code, use the 
  command: make clean; make; make run
 
3_PowerRadiated Contains two folders:
  
  measuredPowerFaradayTensor: -this code calculated the power radiated by
  calculating the electric and magnetic field components of the Faraday tensor for
  an accelerating point charge, as given by the Classical Electrodynamics textbook
  by Jackson. This code was unable to produce adequate results, and it is
  suggested that the theoretical base for the algorithm was not complete in
  considering relativistic phenomena when calculating the Poynting vector.

  measuredPowerLienardWiechert: -this code calculated the power radiated from an
  equation produced in private communications during the project. The equation is
  derived from the Lienard-Wiechert radiation fields of a moving point charge,
  shown in the Classical Electrodynamics textbook by Jackson.

4_LockInAmplifier -the code for the lock-in amplifier can be found here. The
  simulaton setup ran 100 lock-in amplifiers, each scanning for the signal
  frequency within their own bandwidths. 

5_Pictures -pictures used in the latex code for the final report.

In addition to the main codes, each folder contains python/matlab code that was
used to analyse what the issues with code written during the time of the project
were,and to obtain results for the presentation of the work done in the final
report.
