#ifndef DEFS_H    // To make sure you don't declare the function more than once by including the header multiple times.
#define DEFS_H

# include "structure_definitions.h"
const double m =1;
const double inv_mass =1/m;
const double dt=0.01;
const int NrParticles=100;
const double r_cut  = 5;
const double r_cut2 = (r_cut)*(r_cut);
const double sigma =1, epsilon =1; 
const double r_min = pow(2,1/6)*sigma;
const double r_min2= r_min*r_min;
const double rs = 6*r_min/7; // saturation radius, below this potential is assumed linear and force remains constant, to prevent calculation of huge forces at extremely close contacts 
const double rs2=rs*rs;
const double L0 = r_min; // L0 equlibrium length of the dimer molecular spring 
const double Spring_Const=1;
const double red_m = m/2 ; //reduced mass
const double omega = sqrt(Spring_Const/red_m);  // w=sqrt(k/mu) ; spring oscillation frequency 
const double Lx=50, Ly=50, Lz=50;// , R_cut=2.5/*= 1.1225 */,R_shell = 0; // = 0.3775;
const double Volume =Lx*Ly*Lz;
const double Volume_inv = 1/ Volume;
const int cellx=(int) ceil(Lx/r_cut);
const int celly=(int) ceil(Ly/r_cut);
const int cellz=(int) ceil(Lz/r_cut);

const vctr3D box(Lx, Ly, Lz);
const vctr3D rbox(1/Lx, 1/Ly, 1/Lz);
const vctr3D havbox (Lx/2, Ly/2 , Lz/2);
const vctr3D shift(r_cut+r_min,r_cut+r_min,r_cut+r_min);

const double Fs=4*epsilon*(12*pow(sigma/rs,12)-6*pow(sigma/rs,6))/rs;
const double phis =4*epsilon*(pow(sigma/rs,12)-pow(sigma/rs,6));
const double phicut =4*epsilon*(pow(sigma/r_cut,12)-pow(sigma/r_cut,6));
const double sigma6 = sigma*sigma*sigma*sigma*sigma*sigma;
const double sigma12 = sigma6*sigma6;


#endif

