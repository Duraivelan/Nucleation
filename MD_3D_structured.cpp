#include <iostream>
#include <random>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <tuple>
#include <dirent.h>
#include <algorithm>
#include <functional>
#include <array>
# include "defs.h"
# include "force_structured.h"

using namespace std;


// random numbers using random_device option with normal distribution
/*void createInitialPosition_N_particles(std::string fileName, int N) {
      std::random_device rd, rd1;
      std::mt19937 genx(rd()),geny(rd1());
      std::ofstream outFile(fileName);
      std::normal_distribution<> d(0,1);
      for(int i=0;i<N;i++) {
      outFile<<d(genx)<<'\t'<<d(geny)<<std::endl;
      }
      outFile.close();
}
*/

// compare two rows by elements of the first elemtn of the row for sroting 

 struct RowSort {
        bool operator()(vector<int> a, vector<int>  b)
        {   
            return a[0] < b[0];
        }   
    } ;
    
    struct RowUnique {
        bool operator()(vector<int> a, vector<int>  b)
        {   
            return a[1] != b[1];
        }   
    } ;
// random numbers using rand function
void createInitialPosition_N_particles(std::string fileName,std::string fileName2, int N, double Lx, double Ly, double Lz) {
	double x,y,z, vx, vy, vz;
 	srand (time(NULL)); // initialize random seed
 	std::ofstream outFile(fileName);
 	std::ofstream outFile2(fileName2);
 	for(int i=0;i<N;i++) {
 		x=((double) rand() / (RAND_MAX/Lx))-Lx/2;  // create particle position from -Lx/2 to Lx/2
		y=((double) rand() / (RAND_MAX/Ly))-Ly/2;
		z=((double) rand() / (RAND_MAX/Lz))-Lz/2;
		vx= ((double) rand()/(RAND_MAX)-0.5);
		vy= ((double) rand()/(RAND_MAX)-0.5);
		vz= ((double) rand()/(RAND_MAX)-0.5);
		outFile<<x<<'\t'<<y<<'\t'<<z<<std::endl;
		outFile2<<vx<<'\t'<<vy<<'\t'<<vz<<std::endl;
 	}
 	outFile.close();
}
void createInitialPosition_N_dimer_molecules(std::string fileName,std::string fileName2, int N, double Lx, double Ly, double Lz, double L0) {
	double x,y,z, vx, vy, vz, vmolx,vmoly, vmolz;
 	srand (time(NULL)); // initialize random seed
 	std::ofstream outFile(fileName);
 	std::ofstream outFile2(fileName2);
 	for(int i=0;i<N/2;i++) {
 		x=((double) rand() / (RAND_MAX/Lx))-Lx/2;  // create particle position from -Lx/2 to Lx/2
		y=((double) rand() / (RAND_MAX/Ly))-Ly/2;
		z=((double) rand() / (RAND_MAX/Lz))-Lz/2;
		vmolx= ((double) rand()/(RAND_MAX)-0.5);
		vmoly= ((double) rand()/(RAND_MAX)-0.5);
		vmolz= ((double) rand()/(RAND_MAX)-0.5);
		vx= vmolx	+	((double) rand() / (RAND_MAX)) * omega/sqrt(3);
		vy= vmoly	+	((double) rand() / (RAND_MAX)) * omega/sqrt(3);
		vz= vmolz	+	((double) rand() / (RAND_MAX)) * omega/sqrt(3);
		outFile<<x<<'\t'<<y<<'\t'<<z<<std::endl;
		outFile2<<vx<<'\t'<<vy<<'\t'<<vz<<std::endl;
		x=x+L0;  // create particle position at L0 from frist molecule 
		vx=vmolx-(vx-vmolx);
		vy=vmoly-(vy-vmoly);
		vz=vmolz-(vz-vmolz);
		outFile<<x<<'\t'<<y<<'\t'<<z<<std::endl;
		outFile2<<vx<<'\t'<<vy<<'\t'<<vz<<std::endl;
 	}
 	outFile.close();
}
std::vector<int> radialDistFunc(double XYZ[][3], double Lx,double Ly, double Lz, double dr, int N) {
    std::vector<int> rdf((int) floor(sqrt(pow(Lx/2,2)+pow(Ly/2,2)+pow(Lz/2,2)))/dr,0);
	double r;
	for(int j=0;j<N;j++) {
		r=sqrt(pow(XYZ[j][0],2)+pow(XYZ[j][1],2)+pow(XYZ[j][2],2));
	    rdf[(int) floor(r/dr)]+=1;                        // put each particle in a bin according to its position from origin ie. (0,0)
	}
	return rdf;
}

// forceUpdate fucntion included as force.h header file
void verlet( vector<SubData>& particle ) {
	
	for(int i=0;i<NrParticles;i++) 
	{
		particle[i].vel+=particle[i].frc*(0.5*dt*inv_mass);
		particle[i].pos+=particle[i].vel*dt;
		particle[i].pos.PBC(box,rbox);

	}
}

void verletB(vector<SubData>& particle, double vel_scale) {
	if(xxthermo) 
		{
		for(int i=0;i<NrParticles;i++) 
			{
				particle[i].vel+=particle[i].frc*(0.5*dt*inv_mass);
				particle[i].vel=(particle[i].vel)*vel_scale;
			}
       	} 
	else 
		{
		for(int i=0;i<NrParticles;i++) 
			{
				particle[i].vel+=particle[i].frc*(0.5*dt*inv_mass);
			}
       	}
}

int main() {
// current date/time based on current system
   time_t now = time(0);
   struct tm *ltm = localtime(&now);
   cout << "start time"<< '\t'<< ltm->tm_hour << ":";
   cout << ltm->tm_min << ":";
   cout << ltm->tm_sec << endl;
         
int if_create_particles = xxcreate, ifrestart=xxrestart;
         
double kb=1 , T0=0.3, tauT=0.1;
double Temp=0;
double shear_rate = 0; //shear rate
int ifshear = 0;// set equal to 1 for shear
std::string dataFileName="../xxx",dataFileName_new="../xxxnew" ;
int Max_Cluster_N=NrParticles;
double simu_time=dt;
int step=0, nSteps=10000, frame=10;
double vel_scale;
int if_Periodic =1;

std::cout<<cellx<<'\t'<<celly<<'\t'<<cellz<<std::endl;
double  T_Energy, K_Energy, P_Energy, p_energy=0, p_energy_spring=0;
vctr3D L_dimer[NrParticles/2];  // distance between particles of the dimer molecule 
vctr3D dR, dr2;
double R, r2;
double dr=0.05; // step size for RDF calculation
// std::vector<int> RDF((int)  floor(sqrt((Lx/2)*(Lx/2)+(Ly/2)*(Ly/2)+(Lz/2)*(Lz/2)))/dr,0), RDF1((int)  floor(sqrt(Lx*Lx+Ly*Ly))/dr,0);
double KE_rot=0;

vector<SubData>  particle(NrParticles);

// variables for pair detection
const int MaxPairs = 100 ;
const int MaxSplit = 100 ;
int pairs[2][ MaxPairs ][ 3 ] ;		// ! third index: 1,2 = particles, 3 = time of creation / annihilation
int split[2][ MaxSplit ][ 3 ] ;	// ! ibid

int life_span = 10000 	;			// ! life time of pair to qualify as "event"
int save_span = 20000	;			// ! number of steps to be stored per event
int save_step = 0   ; 				// ! frames predating this step are not to be removed
int old_frame;
int ptr_new = 0  ;					// ! pointers, toggle between 1 and 2
int ptr_old = 1  ;

int pairs_old = 0 ;					// ! initiate
int pairs_now = 0 ;
int split_old = 0 ;					// ! initiate
int ancient_birth;
int split_now;
//

if(ifrestart)	{
std::string fileName=dataFileName+"/End_positions.dat";

//read x,y positions from XY.dat
std::ifstream dataFile(fileName);
if(!dataFile.good()) {
	std::cerr<<"Given file is corrupt /n"<<std::endl;
}
else {
	std::cout<<"Reading X, Y, Z co-ordinates"<<std::endl;
    std::string line;
    for (int i=0;i<NrParticles;i++) {
		std::getline(dataFile,line);
    	std::istringstream currentLine(line);    
        currentLine >> particle[i].pos.comp[0];
        currentLine >> particle[i].pos.comp[1];
        currentLine >> particle[i].pos.comp[2];
    }
}	
	fileName=dataFileName+"/Velocities.dat";
	std::ifstream dataFile1(fileName);

	if(!dataFile1.good()) {
	std::cerr<<"Given file is corrupt /n"<<std::endl;
}
else {
	std::cout<<"Reading Vx, Vy, Vz Velocities"<<std::endl;
    std::string line;
     
    for (int i=0;i<NrParticles;i++) {
		std::getline(dataFile1,line);
    	std::istringstream currentLine(line);
        currentLine >> particle[i].vel.comp[0];
        currentLine >> particle[i].vel.comp[1];
        currentLine >> particle[i].vel.comp[2];

    }
}	
} else {

	std::string fileName="../XYZ.dat";
	std::string fileName2="../Initial_Velocities.dat";

if (if_create_particles) {
createInitialPosition_N_dimer_molecules(fileName,fileName2,NrParticles,Lx,Ly,Lz,L0);
}
//read x,y positions from XY.dat
std::ifstream dataFile(fileName);
std::ifstream dataFile2(fileName2);

if(!dataFile.good() && !dataFile2.good() ) {
	std::cerr<<"Given file is corrupt /n"<<std::endl;
}
else {
	std::cout<<"Reading X, Y, Z co-ordinates"<<std::endl;
    std::string line;
    std::string line1;


    for (int i=0;i<NrParticles;i++) {
		std::getline(dataFile,line);
    	std::istringstream currentLine(line);
        currentLine >> particle[i].pos.comp[0];
        currentLine >> particle[i].pos.comp[1];
        currentLine >> particle[i].pos.comp[2];
		std::getline(dataFile2,line1);
    	std::istringstream currentLine1(line1);
        currentLine1 >> particle[i].vel.comp[0];
        currentLine1 >> particle[i].vel.comp[1];
        currentLine1 >> particle[i].vel.comp[2];  
        Temp+=0.5*m*(particle[i].vel.comp[0]*particle[i].vel.comp[0]
				   + particle[i].vel.comp[1]*particle[i].vel.comp[1]
				   + particle[i].vel.comp[2]*particle[i].vel.comp[2]);
              
    }
    	Temp=(Temp)/(1.5*NrParticles*kb);
		vel_scale = sqrt(T0/Temp);
		std::cout<<Temp<<'\t'<<vel_scale<<std::endl;

}

}		
//delete all files before writing data

// following snippet taken from stakcflow link  http://stackoverflow.com/questions/11007494/how-to-delete-all-files-in-a-folder-but-not-delete-the-folder-c-linux
if (ifrestart) {
dataFileName=dataFileName_new;
}
const char *dataFileNamePointer = dataFileName.c_str();  // covnert the datafilename to a char pointer ans pass it to the snippet below which delete all files in that folder before running the simulation
if (!ifrestart) {
struct dirent *next_file;
DIR *theFolder;
char filepath[256];
theFolder = opendir(dataFileNamePointer);
while (( next_file = readdir(theFolder)) )
	{
    // build the full path for each file in the folder
    sprintf(filepath, "%s/%s",dataFileNamePointer, next_file->d_name);
    remove(filepath);
    }
//
}

/* sort particles into cells */
	if (!ifrestart) {

for (int i=0;i<NrParticles;i++) {

	//	particle[i].vel={((double) rand()/(RAND_MAX)-0.5),((double) rand()/(RAND_MAX)-0.5),((double) rand()/(RAND_MAX)-0.5)};
					particle[i].vel=(particle[i].vel)*vel_scale;
	
	
	}
}
for ( int i = 0 ; i < NrParticles/2 ; i ++ ) // 2, not 2.0 , since we want integer value
	{ 
		dR=particle[2*i].pos-particle[2*i+1].pos;
		dR.PBC(box,rbox);
		L_dimer[i]=dR;
	}

std::ofstream outFile(dataFileName+"/K_energy.dat");
std::ofstream outFile1(dataFileName+"/PE_energy.dat");
std::ofstream outFile2(dataFileName+"/Velocities.dat");
std::ofstream outFile7(dataFileName+"/End_positions.dat");
std::ofstream outFile3(dataFileName+"/Pressure_Tensor.dat");
std::ofstream outFile10(dataFileName+"/pairs.dat");

outFile3<<"simu_time"<<'\t'<<"Pxx"<<'\t'<<"Pyy"<<'\t'<<"Pzz"<<'\t'<<"Pxy"<<'\t'<<"Pxz"<<'\t'<<"Pyz"<<'\t'<<"Mxy"<<'\t'<<"Mxz"<<'\t'<<"Myz"<<'\t'<<"Temp"<<std::endl;
// perfrom MD steps
/*	if (ifrestart) {
	simu_time =10.001;
	new_neighbor = TOPMAP(cellx,celly,cellz,if_Periodic,neighbor,box,(R_cut+R_shell),simu_time*shear_rate*Ly);
	std::tie(Pxx[0],Pyy[0],Pzz[0],Pxy[0],Pxz[0],Pyz[0])= forceUpdate(&p_energy ,XYZ, new_neighbor, cell, cellLength, box, (R_cut+R_shell), N, force, shear_rate, simu_time, ifshear , epsilon, sigma, rs);
}
*/
step = 0;
forceUpdate(step, pairs, &pairs_now, ptr_new,  ptr_old, MaxPairs, particle, &p_energy , &p_energy_spring);

simu_time =dt;
do {
	p_energy=0;	
	p_energy_spring=0;
	verlet( particle )	;
	
	pairs_now = 0 ; // ! initiate

 	forceUpdate(step, pairs, &pairs_now, ptr_new,  ptr_old, MaxPairs, particle, &p_energy , &p_energy_spring);
if (pair_detect) {	
	vector<vector<int>> temp_pair(pairs_now+1,vector<int> (3))	;
	
	for (int pn = 1; pn<=pairs_now ; pn++) 
		{ 
		
			
			for (int j = 0; j< 3 ; j ++) 
				{
					temp_pair[pn][j]=pairs[ptr_new	][pn][j];
				}
	//			std::cout<<temp_pair[pn][0]<<'\t'<<temp_pair[pn][1]<<std::endl;

		}		
	//		std::cout<<"before sort "<<'\t'<<pairs_now<<endl;
		
	sort (temp_pair.begin()+1,temp_pair.end(), RowSort());
	
		//		std::cout<<"after sort "<<'\t'<<pairs_now<<endl;
	for (int pn = 1; pn<=pairs_now ; pn++) 
		{ 
	//						std::cout<<temp_pair[pn][0]<<'\t'<<temp_pair[pn][1]<<std::endl;

		}	
	
	
	if(pairs_now>1) {	
	int count=1;
	do
		{
			int j=1;
		do
			{
			if ((temp_pair[count][0]==temp_pair[count+j][0]) && (temp_pair[count][1]==temp_pair[count+j][1]) )
				{
					temp_pair.erase( temp_pair.begin() + count+ j );
					pairs_now-=1;
					j-=1;
				}
				j+=1;
			}	while (j<=(pairs_now-count));
			count=count+1;			
		} while (count<pairs_now);
	}	
	//		std::cout<<"after unique "<<'\t'<<pairs_now<<endl;
	for (int pn = 1; pn<=pairs_now ; pn++) 
		{ 
			for (int j = 0; j< 3 ; j ++) 
				{
					pairs[ptr_new][pn][j]=temp_pair[pn][j];
				}
		//					std::cout<<temp_pair[pn][0]<<'\t'<<temp_pair[pn][1]<<std::endl;

		}	
	
 //	! for all pairs in this step: did they exist before?
  ancient_birth = step - life_span ;
    for (int pn = 1; pn<=pairs_now ; pn++) { 
		for (int po = 1; po<=pairs_old ; po++) { 
			if ( (pairs[	ptr_new	]	[	pn	]	[ 0 ]	== 	pairs[	ptr_old	]	[	po	]	[ 0 ] ) &&
				 (pairs[	ptr_new	]	[	pn	]	[ 1 ] 	== 	pairs[	ptr_old	]	[	po	]	[ 1 ] ) ) 
				{
					pairs[	ptr_new	]	[	pn	]	[ 2 ] = pairs[	ptr_old	]	[	po	]	[ 2 ];	// ! correct time of formation
					pairs[	ptr_old	]	[	po	]	[ 2 ] = step;								// ! mark as processed 
					if ( pairs[	ptr_new	]	[	pn	]	[ 2 ] == ancient_birth ) 
					{		
						save_step = step; 
						outFile10<<pairs[ptr_old][pn][0]<<'\t'<<pairs[ptr_old][pn][1]<<'\t'<<pairs[ptr_old][pn][2]<<'\t'<<"creation"<<std::endl;
						std::cout<<pairs[ptr_old][pn][0]<<'\t'<<pairs[ptr_old][pn][1]<<'\t'<<pairs[ptr_old][pn][2]<<'\t'<<"creation"<<std::endl;

					}
				po=pairs_old;
			}
		} // ! po 
	} // ! pn

// ! for all split pairs: are they still split?
	split_now = 0 ;
	for (int so = 1; so <= split_old; so++) 
	{
		for(int pn = 1; pn <= pairs_now; pn++) 
		{
			if ( ( split[ ptr_old ]	[ so ] [ 0 ] != pairs[	ptr_new	]	[	pn	]	[ 0 ] ) &&
				(  split[ ptr_old ]	[ so ] [ 1 ] != pairs[	ptr_new	]	[	pn	]	[ 1 ] ) ) 
				{
			// ! pair is together again; do not copy to new list
			}
			else
			{
				if ( split[ ptr_old ][ so ] [ 2 ]== ancient_birth ) // ! split for a long time
					{
						save_step = step ;
					   outFile10<<pairs[ptr_new][pn][0]<<'\t'<<pairs[ptr_new][pn][1]<<'\t'<<pairs[ptr_new][pn][2]<<'\t'<<"creation from split"<<std::endl;
					   std::cout<<pairs[ptr_old][pn][0]<<'\t'<<pairs[ptr_old][pn][1]<<'\t'<<pairs[ptr_old][pn][2]<<'\t'<<"creation"<<std::endl;

					}	
				else 
					{
					split_now = split_now + 1 ;
					if ( split_now >= MaxSplit )
					{	
						std::cout<<"Too many splits"<<std::endl;
						abort(); 
					}
				split[ ptr_new ]	[ split_now ] [ 0 ] = split[ ptr_old ]	[ so ] [ 0 ] ;
				split[ ptr_new ]	[ split_now ] [ 1 ] = split[ ptr_old ]	[ so ] [ 1 ] ;
				split[ ptr_new ]	[ split_now ] [ 2 ] = split[ ptr_old ]	[ so ] [ 2 ] ;
				}  // ! ancient / store
			} // ! match
		} // ! pn
	}  // ! so

// ! for all pairs from the previous step: do they qualify as split?
for (int po = 1; po<=pairs_old ; po++) { 
	if ( pairs[	ptr_old	]	[	po	]	[ 2 ] <= ancient_birth ) // then ! old enough to qualify as pair
		{
			split_now += 1 ;
			if ( split_now >= MaxSplit ) 
			{
			    std::cout<<"Too many splits"<<std::endl;
				abort(); // 'Too many splits'
			}
			split[ ptr_new][ split_now][ 0 ] = pairs[ ptr_old][ po][ 0 ];
			split[ ptr_new][ split_now][ 1 ] = pairs[ ptr_old][ po][ 1 ];
			split[ ptr_new][ split_now][ 2 ] = step ; // ! time of break-up
		} //  ! .ne.0
	} // ! po

  pairs_old = pairs_now;
  split_old = split_now;
  ptr_old   = ptr_new;
  ptr_new   = 1 - ptr_new ; // ! toggle
}

	K_Energy=0;
	for ( int i = 0 ; i < NrParticles; i ++ )
			{
						K_Energy+=0.5*m*(particle[i].vel.comp[0]*particle[i].vel.comp[0]
									   + particle[i].vel.comp[1]*particle[i].vel.comp[1]
									   + particle[i].vel.comp[2]*particle[i].vel.comp[2])
									    ;
	}
	
	Temp=(K_Energy)/(1.5*NrParticles*kb);
/*	K_Energy=0;
	for ( int i = 0 ; i < NrParticles/2; i ++ )
			{
				dR=particle[2*i].pos-particle[2*i+1].pos;
				dR.PBC(box,rbox);
				K_Energy=2*0.5*m*((particle[2*i].vel.comp[0]+particle[2*i+1].vel.comp[0])*(particle[2*i].vel.comp[0]+particle[2*i+1].vel.comp[0])*0.25
							+	(particle[2*i].vel.comp[1]+particle[2*i+1].vel.comp[1])*(particle[2*i].vel.comp[1]+particle[2*i+1].vel.comp[1])*0.25
							+	(particle[2*i].vel.comp[2]+particle[2*i+1].vel.comp[2])*(particle[2*i].vel.comp[2]+particle[2*i+1].vel.comp[2])*0.25)
							+	0.5(particle[2*i].vel.comp[0]*particle[2*i].vel.comp[0];
	
			}	*/
if (step%frame==0) { 
		old_frame = step - save_span;
		if ( old_frame > save_step ) 
		{
			char filepath[256];
			sprintf(filepath, "%s",(dataFileName+"/XYZ"+std::to_string((old_frame)/frame)+".xyz").c_str());
			remove(filepath);			
			sprintf(filepath, "%s",(dataFileName+"/VEL"+std::to_string((old_frame)/frame)+".xyz").c_str());
			remove(filepath);
			sprintf(filepath, "%s",(dataFileName+"/pair"+std::to_string((old_frame)/frame)+".dat").c_str());
			remove(filepath);
		}
	
        std::ofstream outFile5(dataFileName+"/XYZ"+ std::to_string(step/frame) +".xyz");
        std::ofstream outFile9(dataFileName+"/VEL"+ std::to_string(step/frame) +".xyz");
        
		outFile5<<NrParticles<<std::endl;
		outFile5<<"X Y Z co-ordinates"<<std::endl;
		// save position, Kinetic energy, Potential energy, Forces every 'frame' steps
		for ( int i = 0 ; i < NrParticles; i ++ )
			{
						outFile5<<'H'<<'\t'<<particle[i].pos.comp[0]<<'\t'<<particle[i].pos.comp[1]<<'\t'<<particle[i].pos.comp[2]<<std::endl;
						outFile9<<particle[i].vel.comp[0]<<'\t'<<particle[i].vel.comp[1]<<'\t'<<particle[i].vel.comp[2]<<std::endl;

	}
      	outFile5<<'\n'<<std::endl;
      	outFile9<<'\n'<<std::endl;
		outFile<<K_Energy<<std::endl;
		outFile1<<p_energy<<'\t'<<p_energy_spring<<std::endl;
		outFile5.close();
		outFile9.close();
	}


	outFile3<<Temp<<'\t'<<KE_rot<<std::endl;
	
	step+=1;
	vel_scale = sqrt(1+(T0/Temp-1)*(dt/tauT));

	verletB( particle , vel_scale) ;

} while(xxnstep);
	

for (int i=0;i<NrParticles;i++) {
	outFile2<<particle[i].vel.comp[0]<<'\t'<<particle[i].vel.comp[1]<<'\t'<<particle[i].vel.comp[2]<<std::endl;
	outFile7<<particle[i].pos.comp[0]<<'\t'<<particle[i].pos.comp[1]<<'\t'<<particle[i].pos.comp[2]<<std::endl;
}
outFile2.close();
outFile.close();
outFile1.close();
outFile3.close();
outFile7.close();

std::ofstream outFile8(dataFileName+"/logfile");
	outFile8<<"NrParticles"<<'\t'<<NrParticles<<std::endl;
	outFile8<<"mass"<<'\t'<<m<<std::endl;
	outFile8<<"kb"<<'\t'<<kb<<std::endl;
	outFile8<<"T0"<<'\t'<<T0<<std::endl;
	outFile8<<"box"<<'\t'<<box.comp[0]<<'\t'<<box.comp[1]<<'\t'<<box.comp[2]<<std::endl;
	outFile8<<"shear rate"<<'\t'<<shear_rate<<std::endl;
	outFile8<<"R_cut"<<'\t'<<r_cut<<std::endl;
	outFile8<<"rs"<<'\t'<<rs<<std::endl;
	outFile8<<"epsilon"<<'\t'<<epsilon<<std::endl;
	outFile8<<"sigma"<<'\t'<<sigma<<std::endl;
outFile8.close();


     // get time now
	now = time(0);
	ltm = localtime(&now);
	cout << "end time"<< '\t'<< ltm->tm_hour << ":";
	cout << ltm->tm_min << ":";
	cout << ltm->tm_sec << endl;
return 0;


}
