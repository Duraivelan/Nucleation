//  Calculate the force between particles i and j

dr=particle[i].pos-(particle[j].pos+dR);
	dr=dr-shift*xxshift;
r2=dr.norm2();

		if (r2<(r_cut2)) 
		{		
			if(pair_detect) {
			*pairs_now += 1;
			if ( *pairs_now >= MaxPairs ) 
			{ 
				step=-1;
				std::cout<<"Too many pairs"<<std::endl;
				abort();
		//	*pairs_now = MaxPairs;
			}

		pairs[	ptr_new	]	[	*pairs_now	]	[ 0 ] = min(i/2,j/2);
		pairs[	ptr_new	]	[	*pairs_now	]	[ 1 ] = max(i/2,j/2);
		pairs[	ptr_new	]	[	*pairs_now	]	[ 2 ] = step;  // ! time of formation
		}
			 if (r2<(rs2)) {
						    	r=std::sqrt(r2);								
												
	// 							exponential potential from PHYSICAL REVIEW E VOLUME 50, NUMBER 3 SEPTEMBER 1994 Browniian dynamics simulations of self-difFusion and shear viscosity of near-hard-sphere colloids
	//							F = n*epsilon*pow(sigma*sigma/rs,n)/rs; //  n is the exponent of the potential function
            					Fij=dr*(Fs/r);
								particle[i].frc+=Fij;
								particle[j].frc-=Fij;

							//	if(cell[neighborList[i][j][k][m][0]][neighborList[i][j][k][m][1]][neighborList[i][j][k][m][2]][lC2] > cell[i][j][k][lC1]) {
								*p_energy+= phis - phicut + Fs*(rs-r);
							} 
			else {				
									
				r2inv=1/r2;
        		r6inv 	= 	r2inv*r2inv*r2inv;
        		r12inv 	= 	r6inv*r6inv;
	//			simple potential
	//			F=2*epsilon*(sigma-r)*rx/r;
	//			*p_energy+=2*epsilon*(sigma-r)*(sigma-r);
	
	//			exponential potential from PHYSICAL REVIEW E VOLUME 50, NUMBER 3 SEPTEMBER 1994 Browniian dynamics simulations of self-difFusion and shear viscosity of near-hard-sphere colloids
	//			F = n*epsilon*pow(sigma*sigma/r2,n/2)/r2; //  n is the exponent of the potential function
	//			F = 4*epsilon*(12*pow(sigma*sigma/r2,6)-6*pow(sigma*sigma/r2,3))/r2; 	WCA potential
				F = 4*epsilon*(12*sigma12*r12inv-6*sigma6*r6inv)*r2inv;
        		Fij=dr*F;
								particle[i].frc+=Fij;
								particle[j].frc-=Fij;

				*p_energy+=4*epsilon*(sigma12*r12inv-sigma6*r6inv) - phicut;

			} 
			

}
