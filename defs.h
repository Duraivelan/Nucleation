#ifndef DEFS_H    // To make sure you don't declare the function more than once by including the header multiple times.
#define DEFS_H

# include "structure_definitions.h"
const double m =1;
const double inv_mass =1/m;
const double dt=0.001;
const double dt2= dt*dt;
const int NrParticles=200;
const double r_cut  = 5;
const double r_cut2 = (r_cut)*(r_cut);
const double sigma =1, epsilon =1; 
const double r_min = pow(2,1/6)*sigma;
const double r_min2= r_min*r_min;
const double rs = 6*r_min/7; // saturation radius, below this potential is assumed linear and force remains constant, to prevent calculation of huge forces at extremely close contacts 
const double rs2=rs*rs;
const double L0 = sigma/5.0; // L0 equlibrium length of the dimer molecular spring 
const double L02 = L0*L0;
const double Spring_Const=10;
const double RTOL = 1.0E-08*L0;
const double RTOL2 = RTOL*RTOL;
const double red_m = m/2 ; //reduced mass
const double omega = sqrt(Spring_Const/red_m);  // w=sqrt(k/mu) ; spring oscillation frequency 
const double Lx=90, Ly=90, Lz=90;// , R_cut=2.5/*= 1.1225 */,R_shell = 0; // = 0.3775;
const double Volume =Lx*Ly*Lz;
const double Volume_inv = 1/ Volume;
const int cellx=(int) ceil(Lx/r_cut);
const int celly=(int) ceil(Ly/r_cut);
const int cellz=(int) ceil(Lz/r_cut);
const double mu = 0.5 ; // reduced mass for a dimer of two identical spheres
const double dt2_by_mu = dt*dt/mu;
const double dt4_by_mu2 = dt2_by_mu*dt2_by_mu;
const vctr3D unit_vec(1.0,1.0,1.0);
const vctr3D box(Lx, Ly, Lz);
const vctr3D rbox(1/Lx, 1/Ly, 1/Lz);
const vctr3D havbox (Lx/2, Ly/2 , Lz/2);
const vctr3D shift(r_cut+r_min,r_cut+r_min,r_cut+r_min);
const int life_span = 50000 ;			// ! life time of pair to qualify as "event"
const int save_span = 100000 ;			// ! number of steps to be stored per event
const int frame=100;
const double Fs=4*epsilon*(12*pow(sigma/rs,12)-6*pow(sigma/rs,6))/rs;
const double phis =4*epsilon*(pow(sigma/rs,12)-pow(sigma/rs,6));
const double phicut =4*epsilon*(pow(sigma/r_cut,12)-pow(sigma/r_cut,6));
const double sigma6 = sigma*sigma*sigma*sigma*sigma*sigma;
const double sigma12 = sigma6*sigma6;


#endif

